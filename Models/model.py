import pandas as pd
import joblib
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, matthews_corrcoef
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier



# Read the data
data = pd.read_csv("dataTain.csv", index_col=0)
print(data.head())

# Assume 'group' is the target variable and the rest are features
X = data.drop(columns=['group'])
y = data['group']

# Encode target variable if necessary (assuming binary classification for example)
y = y.map({'C1': 1, 'C2': 0})

# Define the classifiers to be evaluated
classifiers = {
    "Random Forest": RandomForestClassifier(random_state=42),
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
    "Naive Bayes": GaussianNB(),
    "k-NN": KNeighborsClassifier(),
    "Decision Tree": DecisionTreeClassifier(random_state=42)
}


n_repeats = 5
n_splits = 5
results = []

best_lr_auc = 0
best_lr_model = None

for name, clf in classifiers.items():
    for repeat in range(n_repeats):
        print(f'{name} - Repeat {repeat + 1}/{n_repeats}')
        kfold = KFold(n_splits=n_splits, shuffle=True, random_state=repeat)
        for fold, (train_idx, test_idx) in enumerate(kfold.split(X)):
            print(f' Fold {fold + 1}/{n_splits}')

            # Split the data
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

            # Train the model
            clf.fit(X_train, y_train)

            # Evaluate on the test set
            y_test_pred_prob = clf.predict_proba(X_test)[:, 1]
            y_test_pred = (y_test_pred_prob > 0.5).astype(int)
            y_test_true = y_test

            # Calculate test performance metrics
            test_accuracy = accuracy_score(y_test_true, y_test_pred)
            test_precision = precision_score(y_test_true, y_test_pred, average='binary', zero_division=1)
            test_recall = recall_score(y_test_true, y_test_pred, average='binary')
            test_f1 = f1_score(y_test_true, y_test_pred, average='binary')
            test_auc = roc_auc_score(y_test_true, y_test_pred_prob)
            test_mcc = matthews_corrcoef(y_test_true, y_test_pred)

            print(f'  {name} - Test AUC: {test_auc:.3f}, MCC: {test_mcc:.3f}')

            # Append results to the list
            results.append({
                'Model': name,
                'Repeat': repeat + 1,
                'Fold': fold + 1,
                'Test Accuracy': test_accuracy,
                'Test Precision': test_precision,
                'Test Recall': test_recall,
                'Test F1-Score': test_f1,
                'Test AUC': test_auc,
                'Test MCC': test_mcc
            })

            # Track the best Logistic Regression model
            if name == "Logistic Regression" and test_auc > best_lr_auc:
                best_lr_auc = test_auc
                best_lr_model = clf

# Convert results to DataFrame and display
results_df = pd.DataFrame(results)
print(results_df)

# Save the results to a CSV file
results_df.to_csv('model_evaluation_results.csv', index=False)

# Print the path where the file is saved
print("The model evaluation results have been saved to 'model_evaluation_results.csv'.")

# Save the best Logistic Regression model
if best_lr_model is not None:
    joblib_file = "best_logistic_regression_model.pkl"
    joblib.dump(best_lr_model, joblib_file)
    print(f"The best Logistic Regression model has been saved to '{joblib_file}' with AUC: {best_lr_auc:.3f}.")
