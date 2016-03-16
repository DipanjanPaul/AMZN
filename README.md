Statistical Models for Anomaly Detection

Anomaly Detection is a painful process. Traditional methods involve application of extensive human resource and a cost prohibitive process. Processes get outdated and are seldom updated. Organizations, due to the lack of better choice, have to settle for sub-optimal methods and assessments. This impacts the bottom-line.

Appropriate application of statistical methods and machine learning concepts can improve and optimize anomaly detection. It is also possible to assess the entire population rather than a small sample.

This is an implementation of a customized statistical modeling to quantify the accuracy in a billing process specific to the application of custom discounts and incentives. The incentives are designed for each customer and are customized to provide fixed and variable incentives based on a multiple combinations of service usage. Basically, the incentives are designed for a particular customer for each of their service categories.

The incentives are active for a particular business micro-segment and are predictable (in general), as long as the micro-segments are correctly identified. This model learns the incentive patterns for any given set of customers and their service usage and builds incentive profiles using non-explicit machine learning algorithms. These profiles can be used to certify the billing accuracy of a particular week of billing for the same customers. Using this profile, all billing outliers can be identifies with 99% specificity and 99% sensitivity. 

This model has been developed using R and utilizes the Z statistics and KS principles to identify outliers.
