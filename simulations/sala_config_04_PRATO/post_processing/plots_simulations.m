%% =========================================================================
%  PÓS-PROCESSAMENTO DE RESULTADOS - TCC
%
%  Este script integra a etapa de pós-processamento das simulações numéricas
%  desenvolvidas no âmbito do Trabalho de Conclusão de Curso em Engenharia
%  Acústica. O objetivo é importar os resultados armazenados em arquivo .nmvr,
%  referentes ao caso sala_config_04_PRATO, para posterior análise dos dados
%  e elaboração de gráficos e comparações de interesse.
%
%  Repositório: numav
%  Caso analisado: sala_config_04_PRATO
%  Script: plots_simulations.m
%
%  Autor: Pedro Paulo Calmon Paes
%  Instituição: Universidade Federal de Santa Maria (UFSM)
%
%  Dependência externa:
%  - read_nmvr.m (localizado em ../../../utils)
%
%  Entrada principal:
%  - ../outputs/sala_config_4.nmvr
%
%  Observação importante:
%  Este script pressupõe a manutenção da estrutura de diretórios do projeto,
%  uma vez que utiliza caminhos relativos para localizar tanto a função de
%  leitura quanto o arquivo de resultados.
%% LIMPEZA
clear; clc; close all;

%% ========= IMPORTANDO OS DADOS ============

addpath("../../../utils");
data = read_nmvr('../outputs/sala_config_4.nmvr');

%% Conferência básica das dimensões
disp('Dimensões dos dados importados:')
disp(['freq_stp  : ', mat2str(size(data.freq_stp))]);
disp(['ni2coord  : ', mat2str(size(data.ni2coord))]);
disp(['sei_2_ni  : ', mat2str(size(data.sei_2_ni))]);
disp(['vei_2_ni  : ', mat2str(size(data.vei_2_ni))]);
disp(['cpx_pres  : ', mat2str(size(data.cpx_pres))]);

%% Atalhos
f   = data.freq_stp;
xyz = data.ni2coord;
tri = data.sei_2_ni.';
tet = data.vei_2_ni.';
p   = data.cpx_pres;