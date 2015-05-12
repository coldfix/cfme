#! /usr/bin/bash

commit=$(git rev-parse HEAD)
status=$(git status --untracked-files=no --porcelain)

if [[ -z $status ]]; then
    dirty=true
else
    dirty=false
fi

cat <<EOF

namespace git
{
    const char* commit = "$commit";
    bool has_uncommitted_changes = $dirty;
}

EOF

