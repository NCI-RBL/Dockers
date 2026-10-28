#!/usr/bin/env bash
# Convert BAM <-> CRAM with samtools. Direction is chosen per file from its extension.
# Inputs without coordinate-sorted headers are sorted automatically.
#
# Unix cluster setup (Bash): keep this script in a persistent location and make
# it executable, for example: chmod +x /path/to/cram.sh
# Add this alias to ~/.bashrc, replacing the path with the script's full path:
#   alias CRAM='/path/to/cram.sh'
# Load it in the current shell with `source ~/.bashrc` (or log in again). If the
# cluster uses environment modules, load samtools before running this script;
# also provide a reference name with `-r` or set CRAM_REFERENCE to a name or
# indexed FASTA path.
# Aliases are usually unavailable in scheduler batch scripts; call the script
# by its full path there. The account's home directory should be available on
# both login and compute nodes.
set -uo pipefail

# Defaults can be overridden with the corresponding environment variables.
REFERENCE="${CRAM_REFERENCE:-}"
THREADS="${CRAM_THREADS:-$(nproc 2>/dev/null || echo 4)}"
CRAM_VERSION="${CRAM_VERSION:-3.1}"
OUTDIR=""
FORCE=0
VERIFY=0
DELETE=0
INDEX=1

# Print command syntax, options, and examples.
usage() {
    cat <<EOF
Usage: $(basename "$0") [options] FILE.bam|FILE.cram [...]

Inputs without `SO:coordinate` in the @HD header are coordinate-sorted first.

  *.bam  -> *.cram   (compress)
  *.cram -> *.bam    (decompress)

Options:
    -r REF    Built-in name or FASTA path (indexed with samtools faidx)
                        Built-in names: hg38, mm10                  [env CRAM_REFERENCE: ${REFERENCE:-unset}]
  -t INT    Threads per samtools call                      [env CRAM_THREADS: ${THREADS}]
  -v STR    CRAM version (3.0 or 3.1)                      [env CRAM_VERSION: ${CRAM_VERSION}]
  -o DIR    Output directory                               [default: next to input]
  -f        Overwrite existing outputs
  -n        Do not write an index (.crai/.bai)
  -c        Verify: compare record counts of input and output
  -d        Delete input (and its index) after successful verification (implies -c)
  -h        Show this help

Examples:
    CRAM -r hg38 *.bam
    CRAM -r /path/to/custom.fa sample.cram
  CRAM -t 16 -c sample1.cram sample2.cram
EOF
}

log() { printf '[%s] %s\n' "$(date '+%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# Resolve built-in reference names; add future names to this case statement.
resolve_reference() {
    case "$1" in
        hg38) printf '%s\n' '/fdb/genomebrowser/fasta/hg38/hg38.fa' ;;
        mm10) printf '%s\n' '/fdb/genomebrowser/fasta/mm10/mm10.fa' ;;
        *)    printf '%s\n' "$1" ;;
    esac
}

# Apply command-line options before validating shared prerequisites.
while getopts ":r:t:v:o:fncdh" opt; do
    case "$opt" in
        r) REFERENCE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        v) CRAM_VERSION="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        f) FORCE=1 ;;
        n) INDEX=0 ;;
        c) VERIFY=1 ;;
        d) DELETE=1; VERIFY=1 ;;
        h) usage; exit 0 ;;
        :) die "Option -$OPTARG requires an argument" ;;
        *) usage; die "Unknown option -$OPTARG" ;;
    esac
done
shift $((OPTIND - 1))

[[ $# -gt 0 ]] || { usage; exit 1; }
REFERENCE="$(resolve_reference "$REFERENCE")"
command -v samtools >/dev/null || die "samtools not found in PATH"
[[ -n "$REFERENCE" ]] || die "No reference given (use -r or set CRAM_REFERENCE)"
[[ -s "$REFERENCE" ]] || die "Reference not found: $REFERENCE"
[[ -s "${REFERENCE}.fai" ]] || die "Reference index missing: run 'samtools faidx $REFERENCE'"
[[ "$THREADS" =~ ^[0-9]+$ ]] || die "Threads must be an integer: $THREADS"
[[ -z "$OUTDIR" ]] || mkdir -p "$OUTDIR" || die "Cannot create $OUTDIR"

# Convert each supported input independently and keep a summary of the results.
ok=0; failed=0; skipped=0
failed_files=()

for in in "$@"; do
    case "$in" in
        *.bam)  fmt="cram,version=${CRAM_VERSION}"; ext="cram"; idx="crai" ;;
        *.cram) fmt="bam";                           ext="bam";  idx="bai"  ;;
        *) log "SKIP $in: not .bam or .cram"; skipped=$((skipped + 1)); continue ;;
    esac

    if [[ ! -s "$in" ]] || ! samtools quickcheck "$in"; then
        log "FAIL $in: missing, empty or truncated"
        failed=$((failed + 1)); failed_files+=("$in"); continue
    fi

    base="$(basename "${in%.*}")"
    dir="${OUTDIR:-$(dirname "$in")}"
    out="${dir}/${base}.${ext}"

    if [[ -e "$out" && $FORCE -eq 0 ]]; then
        log "SKIP $in: $out exists (use -f to overwrite)"
        skipped=$((skipped + 1)); continue
    fi

    # Write to a temp name and rename on success so a crash never leaves a valid-looking partial file.
    tmp="${out}.tmp.$$"
    sorted_input="${tmp}.sorted.bam"
    source="$in"
    target="$tmp"
    [[ $INDEX -eq 1 ]] && target="${tmp}##idx##${tmp}.${idx}"
    idx_opt=(); [[ $INDEX -eq 1 ]] && idx_opt=(--write-index)

    log "START $in -> $out"
    if ! header=$(samtools view -H -T "$REFERENCE" "$in"); then
        log "FAIL $in: could not read alignment header"
        failed=$((failed + 1)); failed_files+=("$in"); continue
    fi
    sort_order=$(printf '%s\n' "$header" | awk -F '\t' '$1 == "@HD" { for (i = 2; i <= NF; i++) if ($i ~ /^SO:/) { sub(/^SO:/, "", $i); print $i; exit } }')

    if [[ "$sort_order" != "coordinate" ]]; then
        sort_prefix="${TMPDIR:-$dir}/$(basename "$tmp").sort"
        log "SORT $in: header sort order ${sort_order:-unknown}"
        rm -f "$sorted_input"
        if ! samtools sort -@ "$THREADS" --reference "$REFERENCE" -T "$sort_prefix" \
                -O BAM -o "$sorted_input" "$in" \
            || ! samtools quickcheck "$sorted_input"; then
            log "FAIL $in: coordinate sort failed"
            rm -f "$sorted_input" "$tmp" "${tmp}.${idx}"
            failed=$((failed + 1)); failed_files+=("$in"); continue
        fi
        source="$sorted_input"
    fi

    # Convert to a temporary output, optionally writing the index alongside it.
    if ! samtools view -@ "$THREADS" -T "$REFERENCE" -O "$fmt" ${idx_opt[@]+"${idx_opt[@]}"} \
            -o "$target" "$source" \
        || ! samtools quickcheck "$tmp"; then
        log "FAIL $in: samtools error"
        rm -f "$sorted_input" "$tmp" "${tmp}.${idx}"
        failed=$((failed + 1)); failed_files+=("$in"); continue
    fi
    rm -f "$sorted_input"

    # Optional record-count verification prevents deleting an unverified input.
    if [[ $VERIFY -eq 1 ]]; then
        n_in=$(samtools view -c -@ "$THREADS" -T "$REFERENCE" "$in")
        n_out=$(samtools view -c -@ "$THREADS" -T "$REFERENCE" "$tmp")
        if [[ -z "$n_in" || "$n_in" != "$n_out" ]]; then
            log "FAIL $in: record count mismatch (in=$n_in out=$n_out)"
            rm -f "$tmp" "${tmp}.${idx}"
            failed=$((failed + 1)); failed_files+=("$in"); continue
        fi
        log "VERIFIED $in: $n_in records"
    fi

    mv -f "$tmp" "$out"
    [[ $INDEX -eq 1 ]] && mv -f "${tmp}.${idx}" "${out}.${idx}"
    size_in=$(du -h "$in" | cut -f1)
    size_out=$(du -h "$out" | cut -f1)

    if [[ $DELETE -eq 1 ]]; then
        rm -f "$in" "${in}.bai" "${in}.crai" "${in}.csi" "${in%.*}.bai"
        log "DELETED $in"
    fi

    log "DONE $out ($size_in -> $size_out)"
    ok=$((ok + 1))
done

log "Summary: ${ok} converted, ${skipped} skipped, ${failed} failed"
if [[ $failed -gt 0 ]]; then
    printf '  failed: %s\n' "${failed_files[@]}" >&2
    exit 1
fi
