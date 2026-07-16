#!/bin/bash
GROUP="tga-Kitao-Lab"
SUBMIT_SCRIPT="$(pwd)/submitjob_template.sh"
LOGFILE="submit_11_50_twice_$(date +%Y%m%d_%H%M%S).tsv"


if [ ! -f "$SUBMIT_SCRIPT" ]; then
    echo "错误: 找不到 submitjob_template.sh: $SUBMIT_SCRIPT"
    exit 1
fi

echo -e "trial_id\tjobname1\tjid1\tjobname2\tjid2\tstatus" > "$LOGFILE"

for i in 12; do
    jobname1="Plog_1.$i"
    jobname2="Plog_1.${i}_b"   # 第二段加后缀，避免重名看不清

    echo "=============================="
    echo "提交 TRIAL_ID=$i  ($jobname1 -> $jobname2)"

    # 第一次提交
    out1=$(qsub -g "$GROUP" -N "$jobname1" -v TRIAL_ID="$i" "$SUBMIT_SCRIPT" 2>&1)
    jid1=$(echo "$out1" | grep -oE '[0-9]+' | head -n1)

    if [ -z "$jid1" ]; then
        echo "  第一次提交失败: $out1"
        echo -e "${i}\t${jobname1}\t-\t${jobname2}\t-\tSUBMIT1_FAIL" >> "$LOGFILE"
        continue
    fi
    echo "  第一次提交成功: jid1=$jid1 ($jobname1)"

    # 第二次提交（依赖第一次）
    out2=$(qsub -g "$GROUP" -N "$jobname2" -hold_jid "$jid1" -v TRIAL_ID="$i" "$SUBMIT_SCRIPT" 2>&1)
    jid2=$(echo "$out2" | grep -oE '[0-9]+' | head -n1)

    if [ -z "$jid2" ]; then
        echo "  第二次提交失败: $out2"
        echo -e "${i}\t${jobname1}\t${jid1}\t${jobname2}\t-\tSUBMIT2_FAIL" >> "$LOGFILE"
        continue
    fi
    echo "  第二次提交成功: jid2=$jid2 ($jobname2, hold_jid=$jid1)"

    echo -e "${i}\t${jobname1}\t${jid1}\t${jobname2}\t${jid2}\tOK" >> "$LOGFILE"
done

echo "=============================="
echo "全部提交完成。日志文件: $LOGFILE"
