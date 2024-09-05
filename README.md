# RM-Economterics-QR
# About the project

This project explores different aspects of Quantile Regression as an alternative to traditional Ordinary Least Squares regression, particularly in scenarios where the assumptions of OLS are violated, such as in the presence of heteroscedasticity and outliers. Using Monte Carlo simulations, we assess the quality of the coefficient estimates, obtained using R's quantreg package, over a range of quantiles in terms of finite sample properties. We also assess the quality of predictions using median quantile regression for different settings regarding the underlying data, while comparing it to OLS regression. 

We also address the issue of quantile crossings and provide insight in one mitigation strategy: Log transformation. 

Key findings:
- QR estimates are consistent but slightly biased, with coverage probabilities close to expected confidence levels.
- QR provides more accurate predictions than OLS in datasets with heteroscedasticity and outliers.
- The issue of quantile crossings can be mitigated by applying a log transformation to the data, significantly reducing crossings but not completely eliminating them.

Overall, Quantile Regression proves to be a robust alternative for analyzing complex data, offering a more nuanced view of relationships across the entire data distribution. It is not merely an alternative to OLS but it presents a crucial, nuanced lens through which the multifaceted nature of economic data can be viewed and interpreted. Our findings advocate for a broader adoption, despite its shortcomings, suggesting it is a key to unlocking deeper, more precise economic insights, and guiding more effective policy making.
