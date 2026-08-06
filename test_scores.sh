source ~/.venv/bin/activate

mkdir -p qc

time bin/dummer MET-test.hmm MET-target.fa -S scores.tsv

profiles=$(tail -n +2 scores.tsv | cut -f2 | sort -u)
tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

plot_anchor() {
    local mode="$1" anchor="$2" outbase="$3" title="$4"
    local fams=() files=() prof count i
    for prof in $profiles; do
        count=$(awk -F'\t' -v p="$prof" -v a="$anchor" \
            -v m="$mode" '$1==m && $2==p && $4==a && $5!="-inf" {print $5}' scores.tsv | wc -l)
        if [ "$count" -ge 5 ]; then
            fams+=("$prof")
        else
            echo "# Skipping ${prof} ${mode} ${anchor}: only ${count} finite scores"
        fi
    done
    if [ ${#fams[@]} -eq 0 ]; then
        echo "# Skipping ${mode} ${anchor}: no families with finite scores"
        return
    fi
    local chunk=1
    for ((i = 0; i < ${#fams[@]}; i += 3)); do
        files=()
        for prof in "${fams[@]:i:3}"; do
            awk -F'\t' -v p="$prof" -v a="$anchor" \
                -v m="$mode" '$1==m && $2==p && $4==a {print $5}' scores.tsv \
                > "$tmpdir/${prof}.scores"
            files+=("$tmpdir/${prof}.scores")
        done
        local suffix=""
        [ $chunk -gt 1 ] && suffix="_${chunk}"
        python3 plot_scores.py \
            --title "$title" \
            --fit gumbel \
            --output "qc/${outbase}${suffix}.png" \
            "${files[@]}"
        chunk=$((chunk + 1))
    done
}

for anchor in end start mid; do
    plot_anchor forward-backward "$anchor" "forward_backward_${anchor}" "Forward-backward ${anchor} scores"
done
plot_anchor forward-only all "forward_only_all" "Forward-only scores"
