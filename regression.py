import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

# Load the data
# Replace 'aerodynamic_data.csv' with the path to your file
data = pd.read_csv('Aerodynamics_data-VAI.csv')
data.dropna(inplace=True)

# Extract features (Mach number and angle of attack) and targets (C_l and C_d)
X = data[['Mach', 'AOA']].values
y_cl = data['Lift'].values
y_cd = data['Drag'].values

# Split the data into training and testing sets for both targets
X_train, X_test, y_cl_train, y_cl_test = train_test_split(X, y_cl, test_size=0.2, random_state=42)
_, _, y_cd_train, y_cd_test = train_test_split(X, y_cd, test_size=0.2, random_state=42)

# Define the polynomial degree
degree = 4

# Create polynomial features
poly = PolynomialFeatures(degree=degree, include_bias=False)
X_train_poly = poly.fit_transform(X_train)
X_test_poly = poly.transform(X_test)

# Function to perform polynomial regression and export coefficients
def perform_regression(X_train_poly, y_train, X_test_poly, y_test, target_name):
    # Fit the linear regression model
    model = LinearRegression()
    model.fit(X_train_poly, y_train)

    # Predict on the test set
    y_pred = model.predict(X_test_poly)

    # Evaluate the model
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    # Print the results
    print(f"\nResults for {target_name}:")
    print(f"Mean Squared Error: {mse:.4f}")
    print(f"R^2 Score: {r2:.4f}")

    # Extract feature names and coefficients
    feature_names = poly.get_feature_names_out(['Mach', 'AoA'])
    coefficients = pd.DataFrame({
        'Feature': ['Intercept'] + list(feature_names),
        'Coefficient': np.append(model.intercept_, model.coef_)
    })

    # Save coefficients to a CSV file
    file_name = f'{target_name}_regression_coefficients.csv'
    coefficients.to_csv(file_name, index=False)
    print(f"Coefficients saved to {file_name}")

    return model

# Perform regression for C_l
model_cl = perform_regression(X_train_poly, y_cl_train, X_test_poly, y_cl_test, 'C_l')

# Perform regression for C_d
model_cd = perform_regression(X_train_poly, y_cd_train, X_test_poly, y_cd_test, 'C_d')

mach = 5
aoa = 3
input_data = np.array([[mach, aoa]])
input_data_poly = poly.transform(input_data)

# Predict C_l
cl_prediction = model_cl.predict(input_data_poly)
print(f"Predicted C_l for Mach {mach} and AoA {aoa}: {cl_prediction[0]}")

# Predict C_d
cd_prediction = model_cd.predict(input_data_poly)
print(f"Predicted C_d for Mach {mach} and AoA {aoa}: {cd_prediction[0]}")