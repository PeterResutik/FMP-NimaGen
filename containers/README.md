# Lock files of the Docker image

The image is built from these files, so that every build gets the same tools:

- `FMP-NimaGen.linux-64.lock`: every Conda package of the `FMP-NimaGen` environment, with its version, build and checksum
- `FMP-NimaGen.pip.lock`: the packages pip installs into it (cutadapt, FDSTOOLS and their dependencies)
- `rtn-build.linux-64.lock`: the packages rtn is built with

`FMP-NimaGen.yml` in the main folder names only the main versions, for Conda on any platform.

To change a tool: change `FMP-NimaGen.yml`, build an image from it (in the Dockerfile, create the
environment from `FMP-NimaGen.yml` instead of the lock files), check it, and export the lock files from
that image:

```bash
I=peterresutik/nimagen-pipeline:<version>
docker run --rm --platform linux/amd64 $I micromamba env export -n FMP-NimaGen --explicit --md5 > FMP-NimaGen.linux-64.lock
docker run --rm --platform linux/amd64 $I micromamba env export -n rtn-build --explicit --md5 > rtn-build.linux-64.lock
docker run --rm --platform linux/amd64 $I bash -c 'cd /opt/conda/envs/FMP-NimaGen/lib/python3.12/site-packages &&
  for d in *.dist-info; do [ "$(cat $d/INSTALLER 2>/dev/null)" = pip ] &&
  awk -F": " "/^Name: /{n=\$2} /^Version: /{v=\$2} END{print n \"==\" v}" $d/METADATA; done' | sort -f > FMP-NimaGen.pip.lock
```
