# Motivation

The original brentp/slivar docker container uses BCFtools 1.13, 
which has the plugin split-vep present, but not usable. Attempting 
to run 

`bcftools +split-vep`

causes an error to be thrown, as the bcftools installation is 
unable to use plugins

The split-vep plugin can be crucial in exposing specific 
per-consequence changes to Slivar, which would otherwise require 
custom js code to be written in order to parse these fields.

# Content

This Dockerfile pulls the released slivar container, builds in the
newer bcftools release.. that's all. This exposes the plugins, and
split-vep can be used in the same environment
