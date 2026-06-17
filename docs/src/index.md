# Introduction

We consider a two dimensional optimal control problem in Lagrange form with an affine control system with respect to the (scalar) control. The control problem consists in minimizing the cost functional

```math
\frac{1}{2}\int_{0}^{T} \big(x^2(t) + y^2(t)\big) ~\mathrm dt \to \min,
```
considering that the evolution of the state ``(x,y) \in \mathbb R^2`` is governed by the control system

```math
\displaystyle \dot{x}(t) = u(t), \quad \dot{y}(t) = -u(t) - \frac{y(t)}{2}, \quad |u(t)| \leq M,
```

fixing the initial condition to

```math
\big( x(0), y(0) \big) = (1, -1/2)
```

and considering that the system has to reach the final condition given by

```math
\big( x(T), y(T) \big) = (1/2, 1/2).
```

## Reproducibility

```@setup main
using Pkg
using InteractiveUtils
using Markdown

# Download links for the benchmark environment
function _downloads_toml(DIR)
    link_manifest = joinpath("assets", DIR, "Manifest.toml")
    link_project = joinpath("assets", DIR, "Project.toml")
    return Markdown.parse("""
    You can download the exact environment used to build this documentation:
    - 📦 [Project.toml]($link_project) - Package dependencies
    - 📋 [Manifest.toml]($link_manifest) - Complete dependency tree with versions
    """)
end
```

```@example main
_downloads_toml(".") # hide
```

```@raw html
<details style="margin-bottom: 0.5em; margin-top: 1em;"><summary>ℹ️ Version info</summary>
```

```@example main
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📦 Package status</summary>
```

```@example main
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details style="margin-bottom: 0.5em;"><summary>📚 Complete manifest</summary>
```

```@example main
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```
