#!/bin/bash


cd ../..
pandoc --standalone "doc/documentation/README.md" -c doc/documentation/github.css -o "README.html" --from=gfm  --metadata title="README"
