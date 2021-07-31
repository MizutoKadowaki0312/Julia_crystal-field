# julia-file
.jl repo to save file.


## VSCode で使用するJuliaのversionをアップデートする

1. VSCode > settings を開く．
2. Julia > Executable Path から以下を更新．
```terminal
/Applications/Julia-1.X.app/contents/resources/julia/bin/julia
```





## Julia を Jupyter notebook で使用する

JuliaをアップデートしたらJupyter Notebookで使用するためのコマンド操作．

```julia
pkg > using IJulia
pkg > build IJulia
```
でOK.versionが1.6.0→1.6.1のようなアップデートならこれでOKで，versionが1.6.X→1.7.0のようなアップデートのときは

```julia
pkg > add IJulia
pkg > using IJulia
pkg > build IJulia
```
をする．
