# An-Optimized-deep-learning-framework-for-turbulance-modelling

# 🌀 Turbulent Flow Prediction Using MLP & Cascade MLP Neural Networks

![Python](https://img.shields.io/badge/Python-3.10-3776AB?style=for-the-badge&logo=python&logoColor=white)
![TensorFlow](https://img.shields.io/badge/TensorFlow-2.x-FF6F00?style=for-the-badge&logo=tensorflow&logoColor=white)
![MATLAB](https://img.shields.io/badge/MATLAB-R2022-0076A8?style=for-the-badge&logo=mathworks&logoColor=white)
![Scikit-Learn](https://img.shields.io/badge/Scikit--Learn-1.x-F7931E?style=for-the-badge&logo=scikit-learn&logoColor=white)

> Predicting turbulent eddy viscosity (νt) and Reynolds stress (−u′v′) 
> using MLP and Cascade MLP neural networks trained on DNS data — 
> published as a thesis project.

---

## 📌 Overview

This project integrates **machine learning** with **RANS (Reynolds-Averaged 
Navier–Stokes)** turbulence modeling to improve CFD predictions without 
the cost of full DNS simulations.

Two neural network architectures are implemented:

- **MLP Network** → predicts turbulent eddy viscosity (νt)
- **Cascade MLP Network** → uses MLP output + strain-rate tensor 
  to predict Reynolds stress (−u′v′)

Training is performed on DNS data at **Reτ = 550**, 
and generalization is tested at **Reτ = 180, 390, and 590**.

---

## 🔬 Scientific Background

| Term | Description |
|---|---|
| RANS | Reynolds-Averaged Navier–Stokes equations — industry-standard CFD approach |
| DNS | Direct Numerical Simulation — high-fidelity ground truth data |
| Eddy Viscosity (νt) | Turbulent transport coefficient in RANS models |
| Reynolds Stress (−u′v′) | Turbulent momentum flux — key quantity in wall-bounded flows |
| Reτ | Friction Reynolds number — characterizes flow regime |

---

## 🗂️ Project Structure
turbulence-ml-rans/

│

├── Final Python.ipynb         # MLP & Cascade MLP — training & evaluation

├── Final Matlab.m             # Turbulence data processing & DNS extraction

│

├── requirements.txt           # Python dependencies

│

├── results/

│   └── (plots, prediction vs DNS comparisons)

│

└── README.md

---

## ⚙️ Methodology

### Step 1 — Feature Engineering (MATLAB)
Raw DNS data is processed in MATLAB to extract 10 flow variables:

| Variable | Symbol | Description |
|---|---|---|
| Wall Distance | y | Distance from wall (most important feature) |
| Mean Velocity U | Ū | Streamwise mean velocity |
| Mean Velocity V | V̄ | Wall-normal mean velocity |
| Mean Velocity W | W̄ | Spanwise mean velocity |
| Turbulent Kinetic Energy | TKE | k = 0.5(u′²+v′²+w′²) |
| Strain Rate | S | Rate of deformation tensor |
| Dissipation Rate | ε | Turbulent energy dissipation |
| Cross Correlation | −u′v′ | Reynolds shear stress |
| Eddy Viscosity | νt | Target variable (MLP) |
| Reynolds Stress | τ | Target variable (Cascade MLP) |

### Step 2 — Feature Importance (Python)
A **Random Forest Regressor** ranks input variables by importance.
Wall distance (y) emerged as the most influential feature.
Variables are grouped into input categories and each group is evaluated separately.

### Step 3 — MLP Training (Python)
- Architecture search via **GridSearchCV** (layers: 1–6, neurons: 10–50, optimizers: Adam/SGD/RMSprop)
- Final model: **3 hidden layers** [75 → 65 → 55 neurons], SELU activation
- Optimizer: Adam (lr = 0.001), Loss: MSE
- Training: 20% split, EarlyStopping (patience=16)
- Validation: 5-fold cross-validation

### Step 4 — Cascade MLP (Python)
MLP output (νt) + strain-rate tensor → Cascade MLP → Reynolds stress (−u′v′)

### Step 5 — Generalization Test
Model trained at **Reτ = 550** is tested at **Reτ = 180, 390, 590**
to evaluate cross-Reynolds-number generalization.

---

## 📊 Evaluation Metrics

| Metric | Formula |
|---|---|
| MSE | Mean Squared Error |
| RMSE | Root Mean Squared Error |
| MAE | Mean Absolute Error |
| R² | Coefficient of Determination |

---

## 🛠️ Installation & Usage

### Requirements
```bash
pip install -r requirements.txt
```

`requirements.txt`:numpy, pandas, scipy, tensorflow>=2.10, scikit-learn, matplotlib

### Run
1. Place DNS `.mat` data files in the project root:
   - `data_properties_390.mat`
   - `data_properties_180.mat`
   - `data_properties_590.mat`

2. Open `Final Python.ipynb` in Jupyter Notebook or VS Code

3. Run cells sequentially

---

## 📂 Data

The input data consists of high-fidelity **DNS (Direct Numerical Simulation)** 
results for turbulent channel flow.

- **Training data:** Reτ = 550
- **Test data:** Reτ = 180, 390, 590
- **Format:** `.mat` files (MATLAB format), loaded via `scipy.io.loadmat()`
- **Grid size:** varies by Reynolds number 
  (e.g., 257×256 for Reτ=390, 129×128 for Reτ=180)

> ⚠️ Raw data files are not included due to size constraints.
> DNS data can be obtained from the 
> [Johns Hopkins Turbulence Database](http://turbulence.pha.jhu.edu/) 
> or equivalent sources.

---

## 📈 Key Results

- Wall distance (y⁺) is the most important input feature
- MLP and Cascade MLP both show strong correlation with DNS ground truth
- The model generalizes well across different Reynolds numbers
- Cascade MLP improves Reynolds stress prediction by leveraging νt as intermediate output

---

## 📄 Reference

This project is the implementation of:

> *Prediction of Turbulent Eddy Viscosity and Reynolds Stress Using 
> MLP and Cascade MLP Neural Networks*  
> Master's Thesis -Iran university of science and technology (IUST)

---

## 👤 Behnam Pourpooneh

**[اسمت رو اینجا بذار]**  
M.Sc. in Aerospace engineering  Iran university of science and technology (IUST)
