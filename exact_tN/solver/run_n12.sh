#!/bin/bash
# 一键完成 N=12 决策核验:证明不存在 >69 的族,即 t(12)=69。
# 用法: ./run_n12.sh   (无需参数; 从 0 跑到 4095, 无时间限制)
# 预期输出: "NO clique > 69 with min-index in [0,4095)" => 证毕
# 若输出 "FOUND clique of size 70" => 猜想在 N=12 被推翻(请把输出发给 Claude)
set -e
if [ ! -f ./ck2 ]; then gcc -O2 -o ck2 ck2.c; fi
./ck2 12 69 0 4100 99999999
