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

### **Step 1.1**

Generate the synthetic dataset that will be used to train the neural network:

**Inputs:**

Heston model parameters:

1. delta - volatility of volatility

2. rho - correlation

3. kappa - mean reversion speed
  
4. varsigma - long-term variance

5. v0 - initial variance 

Option's properties:

6. tau - time to maturity

7. S - stock price

8. K - strike price

**Output:**

p - option prices
