# `lib/` — pipeline support classes and dependencies

This directory holds the Groovy helper classes and the one bundled Java dependency that
NovelTree inherits from the [nf-core](https://nf-co.re) pipeline template. Nextflow
automatically adds every `lib/*.groovy` class and every `lib/*.jar` to the classpath at
runtime, so these files are loaded implicitly — none of them is referenced by name in the
pipeline code.

## Contents

| File | Purpose |
| --- | --- |
| `WorkflowMain.groovy` | Startup helpers invoked from `main.nf` (`WorkflowMain.initialize(...)`): help text, parameter-summary logging, config checks, and the pipeline citation string. |
| `NfcoreSchema.groovy` | Validates the launch parameters against [`nextflow_schema.json`](../nextflow_schema.json) and builds the run-parameter summary. Used from `main.nf` (`paramsSummaryMap`) and `WorkflowMain.groovy` (`validateParameters`, `paramsHelp`, `paramsSummaryLog`). |
| `NfcoreTemplate.groovy` | Shared presentation/utility helpers (logo, colours, completion email/summary) used across the two classes above and `main.nf`. |
| `Utils.groovy` | Small utilities (e.g. Conda channel checks); referenced defensively even though Conda is not a supported profile. |
| `nfcore_external_java_deps.jar` | **Load-bearing binary dependency** — see below. |

## `nfcore_external_java_deps.jar` (provenance)

This ~2.2 MB jar is the JSON-schema validation dependency shipped with the nf-core template.
It bundles the Java libraries used by `NfcoreSchema.groovy` to validate parameters — primarily
`org.everit.json.schema` (JSON Schema validation) and `org.json` (JSON parsing), along with
their transitive dependencies (`org.apache.commons.*`, `com.damnhandy.uri.template`,
`com.google.re2j`).

It is **committed intentionally and is required at runtime.** Parameter validation is enabled
by default (`validate_params = true` in `nextflow.config`), and `NfcoreSchema.groovy` imports
these classes directly. Because Nextflow auto-loads `lib/*.jar`, the jar is never referenced by
filename anywhere — but deleting it breaks pipeline startup (the `org.everit`/`org.json` imports
fail to resolve and the whole class fails to compile/load).

Removing the jar would therefore require retiring the entire `NfcoreSchema` validation path
(or migrating it to the modern [`nf-schema`](https://nextflow-io.github.io/nf-schema/) Nextflow
plugin, which keeps `nextflow_schema.json` but replaces this Groovy+jar scaffolding). That is a
deliberate future decision, not a cleanup; for now the jar stays.
