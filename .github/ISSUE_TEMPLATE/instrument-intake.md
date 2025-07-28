---
name: Data Reduction Feature/Request Intake
about: Structured intake form for new feature requests, bugs, or ideas for data reduction software
title: "[INSTRUMENT] Short summary of the request"
labels: instrument-request, needs-triage
assignees: ''
---
Welcome to the Antigen Issue form. Goals of this form are ...
- Help requesters to provide actionable requests to devs
- Encourage completeness without technical overload
- Help fit each request into a collaborative dev workflow
- Automatically tag and milestone issues by instrument
- Automatically notify developers

## Instrument

Which instrument is the request related to?
- [ ] GCMS
- [ ] VIRUSW
- [ ] VIRUS2
- [ ] Other/Describe: `___`

## Type of Request

What kind of request is this? Pick the best match.
- [ ] Bug fix
- [ ] Documentation addition
- [ ] Feature change
- [ ] Feature addition
- [ ] Interface/CLI change
- [ ] Package Changes / Refactor / Dependency / Install
- [ ] Other/Describe: `___`

## Scope of Request

How many things might need to be changed?
- [ ] Small (change to one or a few lines)
- [ ] Medium (change to one file, add a new file)
- [ ] Large (change to many files, breaks interfaces)

## Impact of Request

Help this sort by impact: How large of an effect will this have when completed?
- [ ] Small (one specific use, or minor bug, or current state is usable but not ideal)
- [ ] Medium (general usage, current state is problematic)
- [ ] Large (critical for entire project, current state is very broken or nonexistent)

## Summary

A short summary description of the feature, fix, or idea.
> _Example: "Add support for new CCD channel config used in GCMS for data collection."_

## Motivation and Context

Help provide context: Describe the scientific, operational, or technical motivation.
> _Example: "Without this, data from the upcoming engineering trip cannot be collected."_  

## Constraints or Deadlines

Help prioritize this are related issues:
- [ ] Time-sensitive (e.g. upcoming experiment, upcoming deliverable date)
- [ ] Blocking (e.g. another issue/activity is blocked until this request is completed)
- [ ] Other/Describe: `___`

## Starting Criteria

Things the developer needs to do this:
- [ ] Hardware dependency (e.g. requires use of FCAM, Archon, sensors to test)
- [ ] Data dependency (e.g. requires input data files to test)
- [ ] Host access (e.g. login account)
- [ ] Other/Describe: `___`

Specific data files needed to test:
- [ ] Paths to sample input data (optional): `___`
- [ ] Paths to sample output data as a baseline to compare against (optional): `___`
- [ ] Other/Describe: `___`


## Completion Criteria

Define what a good outcome would look like. Try to frame this as a testable outcome or metric.

> _Example: "Given a folder of input .fits files (paths above), the reduction script writes new FITS output with X-data and new header cards."_

---

## Notifying Watchers

cc @maya-debski @grzeimann @vestuto 


