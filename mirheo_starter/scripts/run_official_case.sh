#!/usr/bin/env bash
set -euo pipefail

ROOT="/home/lzy/projects/mirheo_starter"
source "$ROOT/scripts/activate_mirheo.sh"

name="${1:-}"

case "$name" in
    hello)
        ranks=1
        ;;
    basic|walls|membranes_solvents)
        ranks=2
        ;;
    *)
        echo "用法：bash scripts/run_official_case.sh {hello|basic|walls|membranes_solvents}"
        exit 2
        ;;
esac

src="$MIR_SRC/tests/doc_scripts/${name}.py"
test -f "$src"

mkdir -p "$ROOT/runs"
run_dir="$(mktemp -d "$ROOT/runs/${name}_XXXXXXXX")"

cp "$src" "$run_dir/${name}.py"
cmp "$src" "$run_dir/${name}.py"

if [ "$name" = "membranes_solvents" ]; then
    cp "$MIR_SRC/data/rbc_mesh.off" "$run_dir/rbc_mesh.off"
fi

cd "$run_dir"

sha256sum "$src" "${name}.py" > source_checksums.txt
printf '%s\n' "$run_dir" > run_directory.txt

echo "官方算例：$name"
echo "MPI 进程数：$ranks"
echo "结果目录：$run_dir"

set +e
/usr/bin/time -v \
    /usr/bin/mpirun.openmpi --bind-to none -np "$ranks" \
    "$VIRTUAL_ENV/bin/python" -u "${name}.py" \
    2>&1 | tee console.log
status=${PIPESTATUS[0]}
set -e

printf '%s\n' "$status" > exit_code.txt

if [ "$status" -ne 0 ]; then
    echo "测试未正常结束，退出码：$status"
    echo "请保留日志：$run_dir/console.log"
    exit "$status"
fi

echo "官方脚本正常结束，退出码为 0。"
echo "这表示部署烟雾测试完成，不代表已完成物理精度验证。"
echo "结果目录：$run_dir"
