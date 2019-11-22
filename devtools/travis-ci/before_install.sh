# Temporarily change directory to $HOME to install software
pushd .
cd $HOME
# Make sure some level of pip is installed
python -m ensurepip

if [ "$TRAVIS_OS_NAME" == "osx" ]; then
    HOMEBREW_NO_AUTO_UPDATE=1 brew upgrade pyenv
    # Pyenv requires minor revision, get the latest
    PYENV_VERSION=$(pyenv install --list |grep $PYTHON_VER | sed -n "s/^[ \t]*\(${PYTHON_VER}\.*[0-9]*\).*/\1/p" | tail -n 1)
    # Install version
    pyenv install $PYENV_VERSION
    # Use version for this
    pyenv global $PYENV_VERSION
    # Setup up path shims
    eval "$(pyenv init -)"
fi
pip install --upgrade pip setuptools

# Restore original directory
popd
