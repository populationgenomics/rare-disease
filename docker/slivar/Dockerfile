ARG SLIVAR=brentp/slivar:v0.2.7
FROM $SLIVAR

# update the path to be the same as the slivar base image
ENV PATH="/nim/bin:/root/.nimble/bin:/opt/rh/devtoolset-2/root/usr/bin:$PATH"

# short dockerfile, from the Slivar base, update bcftools to enable plugins
RUN wget https://github.com/samtools/bcftools/releases/download/1.14/bcftools-1.14.tar.bz2 \
    && tar -C /tmp -xf bcftools-1.14.tar.bz2 \
    && rm -f bcftools-1.14.tar.bz2 \
    && cd /tmp/bcftools-1.14 \
    && make prefix=/usr/local install \
    && cd - \
    && rm -rf /tmp/bcftools-1.14;

# this enables the plugin "bcftools +split-vep"
# could add this final command to ensure build works
# RUN bcftools +split-vep

# probably also want python3 in here, but maybe it's a super single purpose image
