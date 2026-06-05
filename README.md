# convert_to_pdbqt

Converte arquivos de ligantes e receptores para o formato **PDBQT**, compatível com PyRx, AutoDock Vina, GNINA e Smina.

---

## Formatos suportados

| Entrada | Tipo | Backend |
|---|---|---|
| `.sdf` | Ligante | Meeko + RDKit |
| `.mol2` | Ligante | Meeko + RDKit |
| `.pdb` | Ligante ou Receptor | Meeko + RDKit / Open Babel |
| `.cif` (mmCIF/PDBx) | Receptor | Open Babel direto ou gemmi + Open Babel |

A detecção automática de tipo (ligante vs. receptor) é feita com base no conteúdo do arquivo. Arquivos `.cif` são sempre tratados como receptores.

---

## Dependências

### Python
```
pip install meeko rdkit gemmi
```

| Pacote | Para que serve |
|---|---|
| `rdkit` | Leitura de moléculas, cargas de Gasteiger, tipos de átomo AutoDock |
| `meeko` | Preparação de ligantes (ligações rotacionáveis, formato ADT) |
| `gemmi` | Fallback para conversão de receptores CIF quando o Open Babel falha |

### Open Babel (externo, obrigatório para receptores)

**Windows:** https://openbabel.org/wiki/Install — marcar "Add to PATH" durante a instalação

**Linux:**
```bash
sudo apt install openbabel
```

**Mac:**
```bash
brew install open-babel
```

---

## Uso

```
python convert_to_pdbqt.py <input> [output] [opções]
```

### Argumentos posicionais

| Argumento | Descrição |
|---|---|
| `input` | Pasta com os arquivos a converter |
| `output` | Pasta de saída (padrão: mesma pasta do input) |

### Opções

| Opção | Valores | Padrão | Descrição |
|---|---|---|---|
| `--type` | `ligand`, `receptor`, `auto` | `auto` | Tipo de molécula. `auto` detecta pelo conteúdo |
| `--formats` | `pdb sdf mol2 cif` | todos | Filtra quais extensões processar |
| `--pH` | número | `7.4` | pH para protonação |
| `--recursive` | — | desligado | Busca arquivos em subpastas |
| `--verbose` | — | desligado | Exibe mensagens detalhadas de debug |

---

## Exemplos

```bash
# Ligantes — detecção automática (SDF, MOL2, PDB)
python convert_to_pdbqt.py ./ligands

# Ligantes — forçar tipo e filtrar formato
python convert_to_pdbqt.py ./ligands --type ligand --formats sdf mol2

# Receptor PDB ou MOL2
python convert_to_pdbqt.py ./proteins --type receptor

# Receptor mmCIF/PDBx (ex: arquivo baixado do RCSB ou AlphaFold)
python convert_to_pdbqt.py ./proteins --type receptor --formats cif

# Pasta de entrada e saída separadas
python convert_to_pdbqt.py ./input ./output

# Busca recursiva em subpastas
python convert_to_pdbqt.py ./all --recursive

# Modo verbose para depuração
python convert_to_pdbqt.py ./proteins --type receptor --verbose
```

---

## Como funciona

### Ligantes (`.sdf`, `.mol2`, `.pdb`)

1. **RDKit** lê o arquivo e adiciona hidrogênios explícitos com coordenadas 3D.
2. **Meeko** identifica ligações rotacionáveis e prepara a topologia de torsões.
3. **RDKit** recalcula as cargas parciais de Gasteiger.
4. O PDBQT é escrito no formato ADT com tipos de átomo AutoDock corretos (`OA`, `NA`, `A`, `HD`, etc.) e registros `CONECT`.

### Receptores `.pdb` / `.mol2`

O **Open Babel** converte diretamente para PDBQT em modo rígido (`-xr`), adicionando apenas hidrogênios polares (`-xh`) e aplicando protonação pelo pH informado.

### Receptores `.cif` (mmCIF/PDBx)

Duas estratégias são tentadas em ordem:

**1. Open Babel direto** *(preferencial)*
O Open Babel suporta mmCIF nativamente desde a versão 2.4. Não gera arquivo intermediário.
```
obabel -icif receptor.cif -opdbqt -O receptor.pdbqt -xr -xh --ph 7.4 -d
```

**2. gemmi + Open Babel** *(fallback automático)*
Usado quando o Open Babel direto falha (versão antiga ou CIF não-padrão). O gemmi lê o mmCIF, remove os hidrogênios existentes, escreve um PDB intermediário limpo, e então o Open Babel converte esse PDB para PDBQT. O PDB intermediário é apagado ao final.

---

## Comportamento do cache

Arquivos já convertidos são ignorados (`[SKIP]`) se o `.pdbqt` existente for mais recente que o arquivo de entrada. Para forçar a reconversão, apague os `.pdbqt` existentes.

---

## Saída esperada

```
========================================================
  Converter -> PDBQT  (PyRx / ADT compatible)
  Mol. type : receptor
  pH        : 7.4
  Input     : C:\Bioinfo\proteins
  Output    : C:\Bioinfo\proteins
  Files     : 3  .cif: 2  .pdb: 1
========================================================
  [CIF] 5HT1A.cif
  [OK] 5HT1A.cif -> 5HT1A.pdbqt  [842 KB, rigid receptor, direct]
  [CIF] 7E2X.cif
  [OK] 7E2X.cif  -> 7E2X.pdbqt   [1204 KB, rigid receptor, direct]
  [OK] protein.pdb -> protein.pdbqt  [310 KB, rigid receptor]
========================================================
  Result: [OK] 3 converted
========================================================
```

Os prefixos de status são:

| Prefixo | Significado |
|---|---|
| `[OK]` | Conversão bem-sucedida |
| `[SKIP]` | Arquivo já convertido e atualizado |
| `[FAIL]` | Erro na conversão |
| `[CIF]` | Pipeline CIF iniciado |

---

## Solução de problemas

**`obabel` não encontrado**
Verifique se o Open Babel está instalado e no PATH. No Windows, certifique-se de marcar "Add to PATH" durante a instalação ou adicionar manualmente `C:\Program Files\OpenBabel-X.X.X` ao PATH do sistema.

**`[FAIL] No atomic sites found in CIF`**
O arquivo CIF pode estar corrompido ou ser um CIF de pequena molécula sem coordenadas fracionárias. Tente abrir o arquivo no Mercury ou VESTA para verificar a integridade. Use `--verbose` para ver a mensagem de erro detalhada do gemmi.

**Ligante sem conformação 3D**
O RDKit tenta gerar coordenadas 3D automaticamente via ETKDGv3. Se falhar, a molécula pode ter problemas de valência. Verifique o arquivo de entrada no RDKit ou no Avogadro.

**Erros de dependência Python**
```bash
pip install --upgrade meeko rdkit gemmi
```
