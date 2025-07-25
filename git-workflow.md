# Git Workflow

- Summary of recommended git workflow with the bare minimum feature branch use 
- Intent is to avoid merge conflicts and to ease backing out of a dev path you later regret.
- Intent is to have `main` branch always working, so keep work-in-progress on feature branches, and do not work on or push directly to `main`.

## GitHub set-up

- Log-in to your GitHub account: `github.com`
- Send a request to be added as a collaborator to [Maya Debski](https://github.com/maya-debski)
- Confirm you can access the repo `Antigen` page here: `https://github.com/maya-debski/Antigen/`
- Set up GitHub authentication, e.g. create and upload SSH key if needed: https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent

## Git CLI set-up

- create a global `~/.gitconfig` for your laptop, for example ...
```bash
[user]
    name = Maya Debski
    email = maya.debski@utexas.edu

[core]
    editor = nano
```

## Git work-flow

- Clone the repo to your laptop:
    - for SSH keys: `git@github.com:maya-debski/Antigen.git`
    - for GitHub access tokens: `git clone https://github.austin.utexas.edu/mcdo-virus2/virus2-docs.git`
- Overview of Workflow
    - Intent is that `main` should always work, or not be a source of merge conflicts, so don't push to `main` directly.
    - Before making edits, create a new feature branch locally, push the feature branch to GitHub, and then open a PR in GitHub and prepend `WIP` to the PR title.
    - Then start editing code and making commits on that feature branch. Push a new commit to GitHub immediately after creating the commit.
    - Leave the `WIP` on the PR for that branch until you are ready for review
    - Let the repo owner merge the branch into `main` once they have reviewed it and approve.
- Workflow command line example
    - `git fetch --all` to get updates from remote
    - `git status` to check for status of current branch. commit, stash, or otherwise backup edits.
    - `get checkout main`
    - `git status`
    - `git pull` # to update `main`
    - `git branch new_branch_name` # create new branch, branching off ff updated `main`
    - `git checkout new_branch_name` # checkout new branch prior to edits
    - `git push -u origin new_branch_name` # push new branch with -u to set up push/pull with github
        - go to GitHub, click "create new PR" for the new branch, possibly have to look under repo "Branches" tab
    - edit a file e.g. `file.py`
    - `git add file.py`
    - `git commit -m "[some key word] do not say did stuff please!"`
    - `git push` # don't have to specify `git push origin new_branch_name` if you're already on the branch, etc
    - iterate until there's a good stopping point, one that works, then...
- On the GitHub web UI for the repo, add comments to the PR discussion thread
- When ready to merge the PR into `main`, remove the `WIP` and add a comment for the repo owner: `Ready to review and merge`
- Owner of repo will merge on github, and then delete the merged branch on github.
- After it's merged, update your local clone/working copy of the repo ...
    - `git fetch --all`
    - `git checkout main`
    - `git status` # check that local main vs remote main are ready for the pull, to avoid local merge conflicts
    - `git pull`  # fetch and merge the new commits from the remote main to your local main
    - if local working copy/clone gets messed up, just `git clone` the repo again but to a different local path
    - if all looks well, and everything was successfully merged into `main`... 
    - `git fetch --prune`  # deletes any remote tracking branches that have been deleted from github.
    - `git branch -D new_branch_name`  # delete the local tracking branch, since all the changes are on main
