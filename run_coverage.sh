#!/bin/bash

#export PATH="/home/jack/.pyenv/bin:$PATH"
#eval "$(pyenv init -)"
#eval "$(pyenv virtualenv-init -)"

coverage run --source bookofnumbers --omit bookofnumbers/factors.py -m pytest tests/cdnf_tests.py
coverage report
coveralls
