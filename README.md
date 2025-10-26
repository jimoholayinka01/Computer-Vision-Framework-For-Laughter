# Computer-Vision-Framework-For-Laughter
This repository contains the code and supporting data for the study: Vision-Based Analysis of Self-Induced Laughter: Affective Dynamics and Facial Action Modelling.”  The research investigates the temporal and affective characteristics of self-induced laughter using computer vision, Facial Action Coding System (FACS), and valence–arousal modelling.
# Vision-Laughter: A Computer Vision Framework for Analyzing Self-Induced Laughter Dynamics

### 🧠 Exploring Affective Dynamics of Self-Induced Laughter Using Computer Vision and Valence–Arousal Modelling

---

## 📘 Project Overview

This repository contains the code and supporting data for the study:  
**“Vision-Based Analysis of Self-Induced Laughter: Affective Dynamics and Facial Action Modelling.”**

The project investigates the **temporal and affective dynamics of self-induced laughter** using **computer vision**, **Facial Action Coding System (FACS)**, and **valence–arousal modelling**.  
By integrating **OpenFace** and **DeepFace**, the framework quantifies frame-level changes in emotional expression and examines residual affective uplift following laughter.

This work contributes to the field of **Affective Computing** and **Digital Mental Health**, providing a reproducible approach for emotion analysis using vision-based models.

---

## 🧩 Repository Structure

VisionLaughter-AffectiveComputing/
│
├── 1_DeepFace_ValenceArousal.ipynb         # Python notebook for DeepFace-based valence–arousal estimation
├── 2_Statistical_Analysis_and_Clustering.R # R script for statistical analysis and k-means clustering
│
├── 3_Plots/                                # Folder containing generated images and plots
│   ├── AU_Trajectories.png
│   ├── Correlation_Heatmap.png
│   └── Cluster_Profiles.png
│
├── 4_Sample_OpenFace_Data.csv              # Example AU data extracted from OpenFace
│
└── README.md                               # Project documentation (this file)


---

## ⚙️ Method Summary

### 1. Data Collection
- 21 participants (aged 18–65) were recorded performing **60 seconds of self-induced laughter**, framed by 5-second pre- and post-laughter segments.  
- Recordings were conducted under natural indoor lighting using a high-definition camera positioned at eye level.  
> ⚠️ **Note:** Original video files are **not included** for **ethical and confidentiality reasons**. Only processed data (CSV format) are provided for reproducibility.

### 2. Feature Extraction
- **Facial Action Units (AUs)** were extracted using **OpenFace 2.0**, focusing on AU06 (cheek raiser), AU12 (lip-corner puller), and AU25 (lips part).  
- **Valence–Arousal** scores were estimated using **DeepFace** on a frame-by-frame basis to derive continuous affective trajectories.

### 3. Statistical and Clustering Analysis
- Conducted in **R**, using *paired t-tests* and *Wilcoxon signed-rank tests* to evaluate pre–post differences in AUs and valence–arousal.  
- **K-means clustering** identified distinct laughter-expression profiles.  
- Visualisations (AU trajectories, correlation heatmaps, and cluster profiles) were created using **ggplot2** and **Matplotlib**.

---

## 📊 Key Findings

- **AU06 and AU12** significantly increased during laughter, confirming their role as **Duchenne laughter markers**.  
- **Valence** increased during laughter and remained elevated afterward, indicating a **residual emotional uplift**.  
- **Clustering** revealed three laughter-expression profiles, suggesting **inter-individual variability** in laughter behaviour.  
- The pipeline demonstrates how **computer vision** can be used as a **lightweight, non-invasive tool** for studying emotional well-being.

---

## 🧠 Technical Requirements

### Python Environment
- Python ≥ 3.10  
- Required Libraries:

### R Environment
- R ≥ 4.3.1  
- Required Packages:

---

## 🔒 Ethical Considerations

All participants provided **informed consent**, and the study was approved by the relevant **institutional ethics committee**.  
Due to **privacy and ethical restrictions**, raw video data are **not publicly shared**.  
Only anonymized, derived data (e.g., AU intensities and valence–arousal CSV files) are included to maintain transparency while protecting participant identity.

---

## 🧩 Future Work

- Extend the dataset with spontaneous laughter comparisons.  
- Incorporate **multimodal data** (audio, EEG, physiological signals).  
- Develop **real-time affective feedback systems** for digital well-being applications.  
- Apply this pipeline to **clinical and stress-monitoring contexts**.

---

## 📚 Citation

If you use this repository or its methods, please cite:

> O. Jimoh, *Vision-Based Analysis of Self-Induced Laughter: Affective Dynamics and Facial Action Modelling*, University of Bolton, 2025.

---

## 💻 Suggested GitHub Repository Name

**`VisionLaughter-AffectiveComputing`**

**Keywords:**  
`Computer Vision` • `Affective Computing` • `DeepFace` • `OpenFace` • `Laughter Analysis` • `Valence–Arousal` • `Mental Health`

---

