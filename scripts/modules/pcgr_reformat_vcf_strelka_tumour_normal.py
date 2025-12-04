#!/usr/bin/env python3

import argparse
import sys
import pysam


def _safe_int(x, default=0):
    try:
        return int(x)
    except Exception:
        return default


def _calc_af_from_ad(sample):
    try:
        ad = sample.get("AD", None)
        if ad and sum(ad) > 0:
            nonref = sum(ad[1:])
            return nonref / float(sum(ad))
    except Exception:
        pass
    return None


def _calc_af_from_snv_counts(sample, ref_base, prefix=""):
    # SNV counts: AU, CU, GU, TU (tumor) or NAU, NCU, NGU, NTU (normal)
    keys = [f"{prefix}AU", f"{prefix}CU", f"{prefix}GU", f"{prefix}TU"]
    if not all(k in sample for k in keys):
        return None
    try:
        au = _safe_int(sample[keys[0]][0])
        cu = _safe_int(sample[keys[1]][0])
        gu = _safe_int(sample[keys[2]][0])
        tu = _safe_int(sample[keys[3]][0])
        total = au + cu + gu + tu
        if total <= 0:
            return None
        ref_map = {"A": au, "C": cu, "G": gu, "T": tu}
        ref_count = ref_map.get((ref_base or "N").upper(), 0)
        alt_count = total - ref_count
        return alt_count / float(total)
    except Exception:
        return None


def _calc_af_from_indel_counts(sample, prefix=""):
    # Indel counts: TAR/TIR (tumor) or NAR/NIR (normal)
    # Some Strelka builds may expose AR/IR as well; try multiple pairs.
    candidates = [
        (sample.get(f"{prefix}TAR"), sample.get(f"{prefix}TIR")),
        (sample.get(f"{prefix}NAR"), sample.get(f"{prefix}NIR")),
        (sample.get(f"{prefix}AR"), sample.get(f"{prefix}IR")),
    ]
    for a, i in candidates:
        try:
            if a is not None and i is not None:
                a0 = _safe_int(a[0])
                i0 = _safe_int(i[0])
                total = a0 + i0
                if total > 0:
                    return i0 / float(total)
        except Exception:
            continue
    return None


def _get_af(sample, record, is_tumor):
    # 1) AF in FORMAT
    af = None
    try:
        if "AF" in sample and sample["AF"]:
            af = float(sample["AF"][0])
    except Exception:
        af = None

    # 2) AD
    if af is None:
        af = _calc_af_from_ad(sample)

    # 3) Strelka counts
    if af is None:
        ref_base = record.ref if record.ref else "N"
        # Try SNV counts first
        af = _calc_af_from_snv_counts(sample, ref_base, prefix="" if is_tumor else "N")
        # If still None, try indel counts
        if af is None:
            af = _calc_af_from_indel_counts(sample, prefix="" if is_tumor else "N")

    # Final fallback
    if af is None:
        af = 0.0
    return float(af)


def reformat_vcf(input_vcf, output_vcf):
    vcf_in = pysam.VariantFile(input_vcf)

    # Ensure the input header itself has the needed INFO tags so records can set them
    header = vcf_in.header
    try:
        header.add_line("##INFO=<ID=TDP,Number=1,Type=Integer,Description=\"Tumor sample depth\">")
    except Exception:
        pass
    try:
        header.add_line("##INFO=<ID=NDP,Number=1,Type=Integer,Description=\"Normal sample depth\">")
    except Exception:
        pass
    try:
        header.add_line("##INFO=<ID=TAF,Number=1,Type=Float,Description=\"Tumor sample AF\">")
    except Exception:
        pass
    try:
        header.add_line("##INFO=<ID=NAF,Number=1,Type=Float,Description=\"Normal sample AF\">")
    except Exception:
        pass

    vcf_out = pysam.VariantFile(output_vcf, "w", header=header)

    for record in vcf_in:
        # Expect two samples: normal first, tumor second in Strelka outputs
        try:
            normal = record.samples[0]
            tumor = record.samples[1]
        except Exception:
            # If sample order unknown, write record unchanged
            vcf_out.write(record)
            continue

        # Depths
        try:
            tdp = _safe_int(tumor.get("DP", 0))
        except Exception:
            tdp = 0
        try:
            ndp = _safe_int(normal.get("DP", 0))
        except Exception:
            ndp = 0
        record.info["TDP"] = tdp
        record.info["NDP"] = ndp

        # Allele fractions
        taf = _get_af(tumor, record, is_tumor=True)
        naf = _get_af(normal, record, is_tumor=False)
        record.info["TAF"] = round(taf, 6)
        record.info["NAF"] = round(naf, 6)

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


def main():
    parser = argparse.ArgumentParser(description="Reformat Strelka VCF for PCGR (add TDP/NDP/TAF/NAF).")
    parser.add_argument("-I", "--input", required=True, help="Input VCF (gz) path")
    parser.add_argument("-O", "--output", required=True, help="Output VCF (gz) path")
    args = parser.parse_args()

    try:
        reformat_vcf(args.input, args.output)
    except Exception as e:
        sys.stderr.write(f"Error during VCF reformat: {e}\n")
        sys.exit(1)


if __name__ == "__main__":
    main()
