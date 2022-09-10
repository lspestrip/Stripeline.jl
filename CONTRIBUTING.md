# How to contribute to the development

Stripeline is a complex project, and we strive to keep it at the best
quality.  This document provides a few guidelines to ensure that the
development is as much organic as possible and that the developers
cooperate harmoniously.

## Version numbers

To make Stripeline as much reliable as possible, we follow semantic
versioning. Therefore, any significant change to Stripeline must thus
be planned and discussed with the lead developers.

## Branches

You should create a separate branch for each feature you want to
implement in the code, and then create a pull request (see below).

## How to use PRs (pull requests)

Start by creating a new branch under `master` with `git checkout
-b`. In this example, we name the branch `my-awesome-feature` and make
it visible on GitHub:

    git checkout master
    git checkout -b my-awesome-feature
    git push --set-upstream origin my-awesome-feature
    
(this example assumes that you named your remote branch on GitHub
`origin`; this is usually the default.)

At this time, you're ready to code! Be sure to commit *frequently*, so
that other people can check how your work is proceeding. Don't worry
if you are creating many small commits—they can be squashed together
once your PR is ready!

To submit the pull request, fire your browser to
https://github.com/lspestrip/Stripeline.jl and pick your branch from
the drop-down list «Branch: …», then click on the button «Pull
request» on the right.


## Things to do before making a PR

- Check that your branch forked from `master` and not from some other
  random branch
- Update `CHANGELOG.md` with a short description of your new features
  under the name `HEAD`.
  
