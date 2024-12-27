FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt update -y && \
    apt install -y make g++ libssl-dev git && \
    apt clean && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/centaur-ai/PWL.git /PWL

WORKDIR /PWL

RUN mkdir deps

RUN git clone https://github.com/asaparov/core.git deps/core && \
    git clone https://github.com/asaparov/math.git deps/math && \
    git clone https://github.com/asaparov/hdp.git deps/hdp && \
    git clone https://github.com/asaparov/grammar.git deps/grammar

RUN make pwl_reasoner_dbg CPPFLAGS_DBG+="-Ideps"

CMD ["/bin/bash"]
