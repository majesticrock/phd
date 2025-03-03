eval "$(ssh-agent -s)"
ssh-add

pull_repo() {
    local repo_path=$1
    echo "Updating $repo_path"
    cd "$repo_path" || exit
    git pull
    cd - || exit
}

pull_repo "ThesisPhd"
pull_repo "PhdUtility"
pull_repo "PresentationsPhd"
pull_repo "phd_plot_scripts"
pull_repo "cpp/ContinuumSystem"
pull_repo "cpp/FermionCommute"
pull_repo "cpp/FlowCommutators"
pull_repo "cpp/Hubbard"
pull_repo "."