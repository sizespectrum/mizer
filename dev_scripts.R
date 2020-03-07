# A file with code snippets useful while developing mizer
# Not part of the mizer package, so added to .Rbuildignore

# Script for upgrading all packaged params objects
for (file in list.files("data")) {
    path <- paste0("data/", file)
    name <- load(path)
    print(name)
    if (is(get(name), "MizerParams")) {
        assign(name, upgradeParams(get(name)))
        save(list = name, file = path, version = 2)
    }
}
