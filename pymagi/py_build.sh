# build python
pip3 install -r pip/requirements.txt
cmake . && make -j $CPU
python3 -c "import pymagi"
if [[ "$1" != "--skip-tests" ]]; then
  nosetests
fi
