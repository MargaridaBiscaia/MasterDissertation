# Calibrating Option Pricing Models using Neural Networks and Population-Based Optimization Methods

**Abstract:**

In finance, model calibration is a crucial task that ensures that option pricing models accurately reflect market conditions, reducing the risk of decisions based on unreliable information. However, this calibration process is often computationally intensive and time-consuming, especially when dealing with complex models.
To address these challenges, neural networks have emerged as a promising approach for developing more efficient option pricing methods, consequently allowing the utilization of algorithms that enhance the calibration process.
We present a novel implementation and comparative analysis of the performance of two types of neural networks, Feedforward Neural Networks (FNN) and Long Short-Term Memory Networks (LSTM), in solving the Heston model calibration problem. We employ a two-step calibration approach, using neural networks to approximate the pricing function and significantly reduce calibration time. Our numerical experiments demonstrate that LSTM networks, particularly when combined with a variant of the Differential Evolution algorithm, can improve calibration accuracy compared to FNNs.

**Keywords:** Neural Networks, Option pricing models, Model's calibration, Population-based optimization methods.

## **Method:**

We adopt the two-step NN-approach, applying it for the calibration of the Heston model. 
The proposed NN-based framework consists on first learning a model and then calibrate it to the data. In the first step, 
a neural network is trained with synthetic data to approximate the pricing function of a model. 
In the second step, the trained network is used in the calibration process instead of the traditional numerical option pricing methods. 

The repository is organized as follows: 

### **Step 1.1**

Generate the synthetic dataset that will be used to train the neural network. This same method is used to generate the out-of-sample dataset and the synthetic calibration dataset.

### **Step 1.2**

Train the FNN and LSTM networks.

### **Step 2**

Calibrate the Heston option pricing model comparing the perfomances of the FNN and LSTM networks, and the ms Differential Evolution and Particle Swarm Optimization algorithms.

### **Data**

Generated datasets for training, out of sample evaluation and synthetic calibration.
