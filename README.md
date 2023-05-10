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

