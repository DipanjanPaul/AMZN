AMZN
====

Statistical Model for Billing Accuracy

This is a statistical model to identify billing accuracy with respect to the incentives applied to the billing data. The incentives are designed for each customer and are customized to provide fixed and variable incentives based on a particular category of service usage. Basically, the incentives are designed for a particular customer for each of their service categories.

The incetives are active for a data range and are predictable (in general). This model learns the inentive patterns for any given set of customers and their service usage and builds incentive profilesusing non explicit machine learning algorithms. These profiles can be used to certify the billing accuracy of a particular week of billing for the same customers. Using these prifile, all billing outliers cna be identifies with a 99% specificity and a 99% sensitivity. 

This model has been developed using R and utilizes the Z statistics and KS principles to identify outliers.

