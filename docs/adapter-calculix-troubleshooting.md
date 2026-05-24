---
title: Troubleshooting the CalculiX adapter
permalink: adapter-calculix-troubleshooting.html
keywords: adapter, calculix, error
summary: "While working with the CalculiX adapter, you may run onto common issues. This is a collection of what we know could go wrong."
---

This list is definitely not complete. If after reading this, you still have issues, please [ask in the preCICE forum](https://precice.discourse.group/).

## Things to check

* Are you using the same version of CalculiX and of the CalculiX adapter? The adapter installation works by replacing files of the original CalculiX, so they should be compatible.
* Can you manually run the `ccx_preCICE` binary?
  * It should be in your `$PATH`
  * If autocompletion does not work (e.g. `ccx_` + TAB key), then it is probably not in your `$PATH`.
* Our tutorials also require CGX (pre- and post-processor of CalculiX).
  * Is CGX installed?
  * Is OpenGL (required by CGX) installed?
