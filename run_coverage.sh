#!/bin/bash

#export PATH="/home/jack/.pyenv/bin:$PATH"
#eval "$(pyenv init -)"
#eval "$(pyenv virtualenv-init -)"

coverage run --source cdnf --omit cdnf/factors.py -m pytest tests/cdnf_tests.py
coverage report
coveralls
