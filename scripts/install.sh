sudo apt update -y
sudo apt install -y make
sudo apt install -y g++
sudo apt-get install -y libssl-dev

git clone https://github.com/centaur-ai/PWL.git

cd PWL
mkdir deps

git clone https://github.com/asaparov/core.git deps/core
git clone https://github.com/asaparov/math.git deps/math
git clone https://github.com/asaparov/hdp.git deps/hdp
git clone https://github.com/asaparov/grammar.git deps/grammar

make pwl_reasoner_dbg CPPFLAGS_DBG+="-Ideps"
