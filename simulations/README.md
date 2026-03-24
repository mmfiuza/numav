\# Simulações acústicas com Numav



Esta pasta contém diferentes configurações de simulações acústicas

utilizando a biblioteca \*\*Numav Solver\*\*.



Cada subdiretório representa um \*\*caso de simulação independente\*\*.



\---



\# Estrutura das simulações



Cada configuração segue o mesmo padrão de organização:

simulations/

└─ sala\_config\_xx/

├─ inputs/

└─ outputs/





\### inputs



Contém os arquivos necessários para executar a simulação:



\- arquivo `.bdf` contendo a malha

\- arquivo `.cpp` que define o problema físico e numérico



Exemplo:

inputs/

sala\_config\_5.cpp

sala\_CONFIG\_5\_versaoteste.bdf





\### outputs



Diretório onde os resultados da simulação são gravados.



Normalmente os resultados são exportados no formato:

.nmvr



utilizado pelas ferramentas de pós-processamento do Numav.



\---



\# Pré-requisitos



Antes de executar qualquer simulação é necessário compilar a biblioteca

\*\*libnumav\*\*.



A partir da raiz do repositório:

rm -rf build &&
cmake -B build -D CMAKE_BUILD_TYPE=Release &&
cmake --build build --parallel ${nproc}


Isso gera a biblioteca:

build/lib/libnumav.a



\---



\# Compilando uma simulação



Cada simulação é compilada manualmente utilizando `g++`

e linkando contra a biblioteca \*\*libnumav\*\*.



Exemplo (configuração 05):

g++ -o simulations/sala_config_05_PRATO/outputs/simul_config_05 simulations/sala_config_05_PRATO/inputs/sala_config_5.cpp -Iinclude -L./build/lib -L/opt/intel/oneapi/mkl/2025.2/lib -Wl,--start-group -l:libmkl_core.a -l:libmkl_gf_ilp64.a -l:libmkl_gnu_thread.a -Wl,--end-group -lgomp -lpthread -lm -ldl -lnumav -m64 -flto

\---



\# Executando a simulação



Após a compilação:



./simul_config_05



\---



\# Resultados



Os resultados são exportados para:

simulations/<config>/outputs/


no formato binário `.nmvr`.



\---



\# Observações



\- Cada pasta de simulação representa um \*\*caso físico diferente\*\*.

\- Os arquivos `.cpp` definem:

  - propriedades físicas

  - fontes acústicas

  - impedâncias de contorno

  - faixa de frequência da simulação.

