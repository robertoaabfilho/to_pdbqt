# Converter to PDBQT

Script em Python para converter arquivos moleculares nos formatos `.pdb`, `.sdf` e `.mol2` para `.pdbqt`, compatível com ferramentas de docking molecular como **PyRx**, **AutoDock Vina**, **GNINA** e **Smina**.

O código diferencia automaticamente, ou por definição do usuário, entre:

- **Ligantes**: convertidos com RDKit + Meeko;
- **Receptores/proteínas**: convertidos com Open Babel.

---

## Funcionalidades

- Conversão de arquivos `.pdb`, `.sdf` e `.mol2` para `.pdbqt`;
- Detecção automática de ligantes e receptores;
- Conversão de ligantes usando **Meeko** e **RDKit**;
- Conversão de receptores usando **Open Babel**;
- Geração de arquivos compatíveis com PyRx e AutoDock Tools;
- Recalcula cargas de Gasteiger para ligantes;
- Atribui tipos atômicos no padrão AutoDock;
- Permite conversão recursiva em subpastas;
- Permite filtrar formatos específicos;
- Evita reconverter arquivos já processados e atualizados.

---

## Formatos suportados

### Entrada

```text
.pdb
.sdf
.mol2
```

### Saída

```text
.pdbqt
```

---

## Requisitos

### Python

Recomenda-se usar Python 3.8 ou superior.

### Dependências Python

Instale as dependências com:

```bash
pip install meeko rdkit
```

Caso esteja usando ambiente virtual:

```bash
python3 -m venv venv
source venv/bin/activate
pip install meeko rdkit
```

No Windows:

```bash
python -m venv venv
venv\Scripts\activate
pip install meeko rdkit
```

---

## Dependência externa: Open Babel

O Open Babel é necessário para converter receptores/proteínas.

### Linux/Ubuntu

```bash
sudo apt update
sudo apt install openbabel
```

### macOS

```bash
brew install open-babel
```

### Windows

Baixe e instale pelo site oficial:

```text
https://openbabel.org/wiki/Install
```

Durante a instalação, marque a opção para adicionar o Open Babel ao `PATH`.

Para testar se o Open Babel foi instalado corretamente:

```bash
obabel -V
```

---

## Estrutura esperada

Você pode organizar os arquivos assim:

```text
projeto/
│
├── convert_to_pdbqt.py
│
├── ligands/
│   ├── composto1.sdf
│   ├── composto2.mol2
│   └── composto3.pdb
│
└── proteins/
    └── receptor.pdb
```

---

## Como usar

### 1. Converter ligantes automaticamente

Se a pasta contém apenas ligantes:

```bash
python3 convert_to_pdbqt.py ./ligands
```

Ou, no Windows:

```bash
python convert_to_pdbqt.py ./ligands
```

Os arquivos `.pdbqt` serão salvos na mesma pasta de entrada.

---

### 2. Converter ligantes especificando o tipo

```bash
python3 convert_to_pdbqt.py ./ligands --type ligand
```

---

### 3. Converter receptores/proteínas

```bash
python3 convert_to_pdbqt.py ./proteins --type receptor
```

---

### 4. Salvar os arquivos convertidos em outra pasta

```bash
python3 convert_to_pdbqt.py ./ligands ./output
```

Exemplo:

```text
Entrada: ./ligands
Saída:   ./output
```

---

### 5. Converter apenas arquivos SDF

```bash
python3 convert_to_pdbqt.py ./ligands --formats sdf
```

---

### 6. Converter apenas arquivos SDF e MOL2

```bash
python3 convert_to_pdbqt.py ./ligands --formats sdf mol2
```

---

### 7. Buscar arquivos também em subpastas

```bash
python3 convert_to_pdbqt.py ./molecules --recursive
```

---

### 8. Definir pH

O pH padrão é 7.4, mas pode ser alterado:

```bash
python3 convert_to_pdbqt.py ./proteins --type receptor --pH 7.0
```

---

### 9. Modo detalhado

```bash
python3 convert_to_pdbqt.py ./ligands --verbose
```

---

## Exemplos práticos

### Converter ligantes para docking no PyRx

```bash
python3 convert_to_pdbqt.py ./ligands --type ligand
```

### Converter receptor para AutoDock Vina

```bash
python3 convert_to_pdbqt.py ./protein --type receptor
```

### Converter todos os arquivos de uma pasta e subpastas

```bash
python3 convert_to_pdbqt.py ./input ./pdbqt_output --recursive
```

### Converter somente moléculas em SDF

```bash
python3 convert_to_pdbqt.py ./sdf_files ./converted --type ligand --formats sdf
```

---

## Como o script funciona

### Para ligantes

O script utiliza:

- **RDKit** para leitura da molécula;
- Adição explícita de hidrogênios;
- Geração de conformação 3D, caso a molécula não tenha coordenadas;
- Cálculo de cargas de Gasteiger;
- **Meeko** para preparar informações relacionadas à flexibilidade e torções;
- Pós-processamento para gerar um `.pdbqt` compatível com AutoDock Tools/PyRx.

O arquivo final inclui:

- Cabeçalho `REMARK`;
- Coordenadas atômicas;
- Cargas parciais;
- Tipos atômicos AutoDock;
- Registros `CONECT`.

---

### Para receptores

O script usa **Open Babel** para gerar o arquivo `.pdbqt`.

A conversão é feita de modo rígido, adicionando hidrogênios polares e aplicando cargas de Gasteiger.

---

## Detecção automática de ligante ou receptor

Quando o parâmetro `--type` não é informado, o script usa o modo:

```bash
--type auto
```

Nesse modo:

- Arquivos `.sdf` e `.mol2` são tratados como ligantes;
- Arquivos `.pdb` são analisados pelo conteúdo;
- Se a maioria dos resíduos detectados for de aminoácidos, o arquivo é tratado como receptor;
- Caso contrário, é tratado como ligante.

Para evitar erros, recomenda-se usar explicitamente:

```bash
--type ligand
```

ou

```bash
--type receptor
```

quando você já souber o tipo da molécula.

---

## Saída esperada

Durante a execução, o script exibe mensagens como:

```text
========================================================
  Converter -> PDBQT  (PyRx / ADT compatible)
  Mol. type : ligand
  pH        : 7.4
  Input     : /caminho/ligands
  Output    : /caminho/ligands
  Files     : 3  .sdf: 3
========================================================

  [OK] composto1.sdf -> composto1.pdbqt  [25 heavy atoms]
  [OK] composto2.sdf -> composto2.pdbqt  [31 heavy atoms]

========================================================
  Result: [OK] 2 converted
========================================================
```

---

## Arquivos já convertidos

Se o arquivo `.pdbqt` já existir e for mais recente que o arquivo original, o script não converte novamente.

Exemplo:

```text
[SKIP] composto1.sdf (already converted)
```

Isso evita processamento repetido desnecessário.

---

## Possíveis erros e soluções

### Erro: dependência Python ausente

Mensagem possível:

```text
[FAIL] Missing dependency
```

Solução:

```bash
pip install meeko rdkit
```

---

### Erro: Open Babel não encontrado

Mensagem possível:

```text
[FAIL] 'obabel' not found.
```

Solução no Ubuntu:

```bash
sudo apt install openbabel
```

Depois teste:

```bash
obabel -V
```

---

### Erro ao ler molécula

Mensagem possível:

```text
[FAIL] Could not read: arquivo.sdf
```

Possíveis causas:

- Arquivo corrompido;
- Formato inválido;
- Molécula sem estrutura química válida;
- Problemas de sanitização pelo RDKit.

Sugestão:

- Verifique o arquivo original;
- Teste abrir em Avogadro, PyMOL, Discovery Studio ou Open Babel;
- Tente converter o arquivo para outro formato antes de rodar o script.

---

### Erro de timeout

Mensagem possível:

```text
[FAIL] Timeout processing receptor.pdb
```

O script possui limite de 300 segundos para conversão com Open Babel.

Possíveis soluções:

- Verificar se o receptor é muito grande;
- Remover águas, ligantes e cadeias desnecessárias antes da conversão;
- Separar cadeias proteicas;
- Preparar o receptor previamente em outro software.

---

## Recomendações antes do docking

Antes de usar os arquivos `.pdbqt` em docking molecular, recomenda-se:

### Para ligantes

- Conferir se a estrutura 3D está correta;
- Verificar protonação;
- Verificar cargas;
- Conferir se a molécula não possui fragmentos desconectados;
- Confirmar se o arquivo `.pdbqt` abre corretamente no PyRx ou AutoDock Tools.

### Para receptores

- Remover moléculas de água, se necessário;
- Remover ligantes cristalográficos, se não forem usados;
- Corrigir resíduos ausentes;
- Adicionar hidrogênios;
- Conferir protonação em pH fisiológico;
- Verificar se o receptor foi tratado como rígido.

---

## Comando geral

```bash
python3 convert_to_pdbqt.py INPUT_FOLDER [OUTPUT_FOLDER] [opções]
```

### Argumentos principais

| Argumento | Descrição |
|---|---|
| `input` | Pasta contendo arquivos `.pdb`, `.sdf` ou `.mol2` |
| `output` | Pasta de saída. Opcional. Se omitida, salva na pasta de entrada |
| `--type` | Tipo molecular: `ligand`, `receptor` ou `auto` |
| `--formats` | Formatos a converter: `pdb`, `sdf`, `mol2` |
| `--pH` | pH usado na protonação. Padrão: `7.4` |
| `--recursive` | Busca arquivos também em subpastas |
| `--verbose` | Exibe mensagens mais detalhadas |

---

## Exemplo completo

```bash
python3 convert_to_pdbqt.py ./input_molecules ./output_pdbqt --type ligand --formats sdf mol2 --recursive --pH 7.4
```

Esse comando:

- Lê arquivos da pasta `input_molecules`;
- Procura também em subpastas;
- Converte apenas `.sdf` e `.mol2`;
- Trata todos os arquivos como ligantes;
- Salva os `.pdbqt` em `output_pdbqt`;
- Usa pH 7.4.

---

## Aplicação

Este script é útil para preparar arquivos de entrada para estudos de docking molecular em:

- PyRx;
- AutoDock Vina;
- AutoDock Tools;
- GNINA;
- Smina;
- pipelines automatizados de triagem virtual.

---

## Observação importante

Embora o script automatize a conversão para `.pdbqt`, a preparação molecular para docking pode exigir curadoria manual, especialmente para receptores proteicos. Sempre valide visualmente os arquivos finais antes de iniciar o docking.

---

## Autor

Script para conversão automatizada de estruturas moleculares para `.pdbqt`, com foco em compatibilidade com PyRx, AutoDock Tools e ferramentas relacionadas de docking molecular.
