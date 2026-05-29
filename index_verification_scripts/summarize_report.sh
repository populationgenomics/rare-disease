#!/usr/bin/env bash
# summarize_report.sh <report_file>

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 <report_file>"
  exit 1
fi

REPORT_FILE="$1"
SUMMARY_FILE="summary_report.txt"

# 1. Count the total, OK, and TRUNCATED samples
total_count=$(grep -E "^(OK|TRUNCATED)" "$REPORT_FILE" | wc -l)
ok_count=$(grep -c "^OK" "$REPORT_FILE")
truncated_count=$(grep -c "^TRUNCATED" "$REPORT_FILE")

# 2. Write the summary data to the final file
{
  echo "========================================"
  echo "          GENOMICS RUN SUMMARY          "
  echo "========================================"
  echo "Total Samples Processed : $total_count"
  echo "OK Samples (Healthy)    : $ok_count"
  echo "Truncated Samples       : $truncated_count"
  echo "----------------------------------------"
  echo "Truncated File Paths:"
  echo "----------------------------------------"

  # 3. Extract just the file paths for truncated samples (column 2)
  if [ "$truncated_count" -gt 0 ]; then
    grep "^TRUNCATED" "$REPORT_FILE" | awk -F'\t' '{print $2}'
  else
    echo "No truncated samples found! All clear."
  fi
  echo "========================================"
} > "$SUMMARY_FILE"

echo "Summary completed! Results written to: $SUMMARY_FILE"
