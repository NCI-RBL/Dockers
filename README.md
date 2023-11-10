[![DOI](https://zenodo.org/badge/576390595.svg)](https://zenodo.org/doi/10.5281/zenodo.10107622)



# Dockers


A repository of recipes for Docker containers built for the NCI-RBL.


Build and Publish Docker container:

```bash
bash ./build Dockerfile.vX.X wilfriedguiblet/containername:vX.X
docker push wilfriedguiblet/containername:vX.X
```

Build and Publish GitHub Page:

```bash
mkdocs build
mkdocs gh-deploy
```

# Workflow documentation

For documentation on how to run the different worflows, please visit the [GitHub Page](https://nci-rbl.github.io/Dockers/).
