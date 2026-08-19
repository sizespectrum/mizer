# How mizer code tells the user something

Use this skill whenever you are about to add a `message()`, a `warning()`, or
any other report to mizer: a default that was filled in, an input that was
adjusted, an instruction that could not be carried out.

Nearly everything mizer says while setting up or changing a model goes through
the signal mechanism in `R/info_signals.R` — 61 `signal_*()` calls collected at
13 `with_info_level()` sites, against 9 surviving plain `message()` calls, all
of them deliberate. A plain `message()` added to a setter looks correct, passes
its tests, and is wrong: it ignores `info_level`, it is not collected with the
other reports, and on the `species_params<-()` path it is swallowed. That is
what #489 was, and moving the stragglers into the mechanism is most of what PR
#508 did.

The mechanism itself is documented in the roxygen of `with_info_level()` and
`signal_info()`. Read those for how it works. This skill is about which call to
make, and about the two choices that are silently wrong.

## Which call to make

| Situation | Call |
|---|---|
| Any new report about a choice mizer made or an input it adjusted | `signal_info()` |
| A setter is leaving a frozen array alone although the parameters say otherwise | `signal_not_recalculated()` |
| The user changed something and it is being ignored because an array is frozen | `signal_frozen_changes()` — do not raise `signal_frozen()` directly |
| Progress during a long job ("Computing …") | plain `message()` |
| Deprecation | `lifecycle::deprecate_warn()` |
| Invalid argument, `gear_params` row checks | `stop()` |

The last three are deliberately outside the mechanism: collected reports are
given at the *end* of the call, which is no use for progress, and an error has
to stop the work rather than be gathered up. Everything else belongs inside.

New user-facing entry point that reports? Give it
`info_level = default_info_level()` and wrap its body:

```r
myFunction <- function(params, ..., info_level = default_info_level()) {
    with_info_level(info_level = info_level, {
    ...
    })
}
```

Handlers nest by themselves — an inner one steps aside — so wrap without
checking what your caller did.

## The two choices that are silently wrong

**`severity` is not a judgement about how serious the report is.** It is a fact
about whether the report has to survive `suppressMessages()`. `species_params<-()`
runs `suppressMessages()` over its recalculation to quieten the routine chatter,
so a `severity = "info"` report raised anywhere under it never reaches the user.

> Use `severity = "warning"` when the user asked for something that is not
> happening. Use `"info"` when mizer is telling them about a choice it made.

Reasoning from "is this bad enough to warrant a warning?" reproduces #489: the
report exists, the tests pass, nobody sees it.

**`unhandled` decides what happens when nothing is collecting**, for instance
when a rate setter is called directly rather than through `setParams()`. Both
values are invisible at the call site:

- `"drop"` — say nothing. Right for chatter that only makes sense as part of a
  report about a whole model.
- `"show"` — report it there and then. Right when this may be all the user
  hears. Every report about an instruction that is not being carried out uses
  it.

`level` is the mild one: 1 survives `info_level = 1`, 3 is chatter that only the
default shows. When in doubt, an instruction that is not being carried out is
level 1.

## Registering a new freezable rate array

If you add a slot that can be frozen — set by hand and thereafter protected by a
comment — add an entry to `frozen_rate_params()` giving the quantity as the user
knows it, the `reset = TRUE` call that unfreezes it, and the parameters that
feed it. Without an entry the user gets only the setter's message and no
warning.

**A wrong entry is worse than a missing one.** A missing parameter costs a
warning; a parameter listed that does not in fact feed the array warns about a
change that *did* take effect. List what the setter reads directly plus the main
inputs of the default calculations, and no more.

This is the third registration a new species parameter may need — see
`.claude/skills/species-param-defaults.md` for the other two.

Never decide "the user's change was ignored" by comparing the frozen array
against what the formula gives. Mizer freezes `mu_b`, `cc_pp` and `psi` itself
when it builds the trait-based and community models, so that test warns on every
parameter change in those models, including a no-op assignment. Decide it from
the parameters that changed, which is what `signal_frozen_changes()` does.

## Testing

Reports are given **on exit**, collapsed by content, and joined with `"\n"`, so
the idioms differ from testing a plain `message()`:

```r
# A message ends in a newline; a warning does not.
expect_message(with_info_level(signal_info("a", "note", level = 1)), "^note\n$")
expect_warning(with_info_level(signal_frozen("metab", "frozen")), "^frozen$")

# Messages and warnings are reported separately, so nest the expectations.
expect_warning(expect_message(with_info_level({...}), "^info\n$"), "^frozen$")

# Set the option with withr so it is restored.
withr::local_options(mizer_info_level = 0)
```

- `expect_message()` on one report inside a wrapped call **fails** — it arrives
  at the end, joined with everything else raised in the same call. Match on the
  combined text, or test the `signal_*()` function on its own.
- A `severity = "warning"` report needs `expect_warning()`, however routine it
  reads.
- Identical reports collapse to one line; two different things said about the
  same `var` are both kept. Assert on both if a change should produce both.
- Tests for the frozen warning belong in `test-species_params.R` — the file
  named after the R file defining the function under test, not after the
  feature. The mechanism's own tests are in `test-info_signals.R`.

## Documenting it

A new report that can appear on code a user has already written is a behaviour
change: add it to `NEWS.md` and to
`inst/skills/upgrade-mizer-code/references/mizer-<version>.md` for the release
being prepared, with a row in that skill's symptom index, then regenerate with
`build_guides()`. Keep the entry to what an upgrader needs — the design
rationale belongs in `NEWS.md` alone.
