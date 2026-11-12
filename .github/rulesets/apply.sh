#!/usr/bin/env bash
# Apply versioned repository settings and rulesets. Requires gh-axi with admin access.
# Run again after changing a policy; existing rulesets are updated by name.
# The `develop` branch must exist before the default branch can be set.
set -euo pipefail

repository="${1:-matthiaskoenig/porous_media}"
rules_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

command -v gh-axi >/dev/null || { echo 'Install and authenticate gh-axi first.' >&2; exit 1; }

gh-axi api PATCH "repos/${repository}" \
  --field default_branch=develop \
  --field allow_auto_merge=true \
  --field delete_branch_on_merge=true \
  --field allow_update_branch=true \
  --field allow_squash_merge=true \
  --field allow_rebase_merge=true \
  --field allow_merge_commit=true \
  --jq '{default_branch,allow_auto_merge,delete_branch_on_merge,allow_update_branch,allow_squash_merge,allow_rebase_merge,allow_merge_commit}'

for policy_path in "${rules_directory}"/*.json; do
  policy_name="$(basename "${policy_path}" .json)"
  policy_id="$(gh-axi api "repos/${repository}/rulesets" --jq "map(select(.name == \"${policy_name}\")) | .[0].id // empty")"
  if [[ -n "${policy_id}" ]]; then
    gh-axi api PUT "repos/${repository}/rulesets/${policy_id}" --input "${policy_path}" --jq '{id,name,enforcement}'
  else
    gh-axi api POST "repos/${repository}/rulesets" --input "${policy_path}" --jq '{id,name,enforcement}'
  fi
done
