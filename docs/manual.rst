Manual pySTARMA
================
This file contains the manual for using the pySTARMA library

SPACE TIME ARMA (STARMA Object)
-----------------
Description
~~~~~~~~~~~~~~~~~~~~~~
The **STARMA class** can be used to estimate **STARMA models**. The method ``STARMA.fit()`` performs the estimation of the model parameters. The method ``STARMA.predict()`` executes the forecast (still in the development stage). The method ``STARMA.get_model()`` returns the full model. The ``STARMA.get_item()`` method returns a selected property of the model (`see Return Values STARMA`_).

Usage
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python
  
  model = sm.STARMA(p, q, ts_matrix, wa_matrices, iterations(optional))
  model.fit()
  model.get_model()
  model.get_item()
  
Example  
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

  from pySTARMA import starma_model as sm
  
  #Create instance of STARMA
  model = sm.STARMA(5, 2, time_series, wa_matrices, 3)
  
  #Estimate parameters
  model.fit()

  #Print explicit item 
  print(model.get_item('bic'))
  
Attributes
~~~~~~~~~~~~~~~~~~~~~~
+---------------------+---------------------------------------------+
| Attribute           | Value                                       |
+=====================+=============================================+
|p                    | Number or list of autoregressive parameters |
+---------------------+---------------------------------------------+
|q                    | Number or list of moving average parameters |
+---------------------+---------------------------------------------+
|ts_matrix            | Time series matrix                          |
+---------------------+---------------------------------------------+
|wa_matrices          | List of adjacency matrices                  |
+---------------------+---------------------------------------------+
|iterations(optional) | Number of iteration of kalman filtering,    |
|                     | only for estimation of moving average       |
|                     | parameters                                  |
+---------------------+---------------------------------------------+

Return Values
~~~~~~~~~~~~~~~~~~~~~~

.. _`see Return Values STARMA`:

A dictionary is returned as a 'model' with the following values:

+---------------------+---------------------------------------------+
| Value               | Description                                 |
+=====================+=============================================+
|residuals            | Matrix with estimated residuals             |
+---------------------+---------------------------------------------+
|phi                  | Matrix with estimated AR-parameters         |
+---------------------+---------------------------------------------+
|phi_tvalue           | Matrix with estimated AR-t-values           |
+---------------------+---------------------------------------------+
|phi_pvalue           | Matrix with estimated AR-p-values           |
+---------------------+---------------------------------------------+
|theta                | Matrix with estimated MA-parameters         |
+---------------------+---------------------------------------------+
|theta_tvalue         | Matrix with estimated MA-t-values           |
+---------------------+---------------------------------------------+
|theta_pvalue         | Matrix with estimated MA-p-values           |
+---------------------+---------------------------------------------+
|sigma2               | Standard deviation                          |
+---------------------+---------------------------------------------+
|bic                  | Bayesian information criterion              |
+---------------------+---------------------------------------------+

SPACE TIME ARIMA (STARIMA Object)
-----------------

Description
~~~~~~~~~~~
The **STARIMA class** can be used to estimate **STARIMA models**. The method ``STARIMA.fit()`` performs the estimation of the parameters. The method ``STARIMA.predict()`` executes the forecast (still in the development stage). The method ``STARIMA.get_model()`` returns the full model. The ``STARIMA.get_item()`` method returns a selected property of the model (`see Return Values STARIMA`_).

Usage
~~~~~~
.. code-block:: python
  
  model = sm.STARIMA(p, q, d, ts_matrix, wa_matrices, iterations(optional))
  model.fit()
  model.get_model()
  model.get_item()
  
Example  
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python

  from pySTARMA import starma_model as sm
  
  #Create instance of STARIMA
  model = sm.STARMA(5, 2, (1,), time_series, wa_matrices, 3)
  
  #Estimate parameters
  model.fit()

  #Print explicit item 
  print(model.get_item('bic'))
  
Attributes
~~~~~~~~~~~~~~~~~~~~~~
+---------------------+---------------------------------------------+
| Attribute           | Value                                       |
+=====================+=============================================+
|p                    | Number or list of autoregressive parameters |
+---------------------+---------------------------------------------+
|q                    | Number or list of moving average parameters |
+---------------------+---------------------------------------------+
|d                    | List of numbers of differentiations         |
+---------------------+---------------------------------------------+
|ts_matrix            | Time series matrix                          |
+---------------------+---------------------------------------------+
|wa_matrices          | List of adjacency matrices                  |
+---------------------+---------------------------------------------+
|iterations(optional) | Number of iteration of kalman filtering,    |
|                     | only for estimation of moving average       |
|                     | parameters                                  |
+---------------------+---------------------------------------------+

Return Values
~~~~~~~~~~~~~~~~~~~~~~

.. _`see Return Values STARIMA`:

A dictionary is returned as a 'model' with the following values:

+---------------------+---------------------------------------------+
| Value               | Description                                 |
+=====================+=============================================+
|residuals            | Matrix with estimated residuals             |
+---------------------+---------------------------------------------+
|phi                  | Matrix with estimated AR-parameters         |
+---------------------+---------------------------------------------+
|phi_tvalue           | Matrix with estimated AR-t-values           |
+---------------------+---------------------------------------------+
|phi_pvalue           | Matrix with estimated AR-p-values           |
+---------------------+---------------------------------------------+
|theta                | Matrix with estimated MA-parameters         |
+---------------------+---------------------------------------------+
|theta_tvalue         | Matrix with estimated MA-t-values           |
+---------------------+---------------------------------------------+
|theta_pvalue         | Matrix with estimated MA-p-values           |
+---------------------+---------------------------------------------+
|sigma2               | Standard deviation                          |
+---------------------+---------------------------------------------+
|bic                  | Bayesian information criterion              |
+---------------------+---------------------------------------------+



Space Time Autocorrelation Function (STACF Object)
-----------------

Description
~~~~~~~~~~~~~~~~~~~~~~
With the **STACF class**, the space-time-autocorrelation-function can be estimated.

Usage 
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python
 
  stacf = Stacf(ts_matrix, wa_matrices, t_lags)
  stacf.estimate()
  stacf.get()

Example
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python
  
  from pySTARMA import stacf_stpacf as st
  
  #Create instance of STACF
  stacf = st.Stacf(time_series, weight_matrices, 25)

  #Estimate STACF
  stacf.estimate()

  #Print estimated STACF
  print(stacf.get())

Attributes
~~~~~~~~~~~~~~~~~~~~~~
+---------------------+---------------------------------------------+
| Attribute           | Value                                       |
+=====================+=============================================+
|ts_matrix            | Time series matrix                          |
+---------------------+---------------------------------------------+
|wa_matrices          | List of adjecency matrices                  |
+---------------------+---------------------------------------------+
|t_lags               | Number of time lags                         |
+---------------------+---------------------------------------------+

Return Values
~~~~~~~~~~~~~~~~~~~~~~
List with lists for each spatial lag. Spatial lags lists contains the estimated spatial autocorrelation for the corresponding time lag. 
  
  List index 0 --> time lag 0 etc..

Space Time Partial Autocorrelation Function (STPACF-Object)
-----------------

Description
~~~~~~~~~~~~~~~~~~~~~~
With the **STPACF class**, the space-time-partial-autocorrelation-function can be estimated.

Usage 
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python
 
  stpacf = Stpacf(ts_matrix, wa_matrices, t_lags)
  stpacf.estimate()
  stpacf.get()

Example
~~~~~~~~~~~~~~~~~~~~~~
.. code-block:: python
  
  from pySTARMA import stacf_stpacf as st
  
  #Create instance of STACF
  stpacf = st.Stpacf(time_series, weight_matrices, 25)

  #Estimate STACF
  stpacf.estimate()

  #Print estimated STACF
  print(stpacf.get())

Attributes
~~~~~~~~~~~~~~~~~~~~~~
+---------------------+---------------------------------------------+
| Attribute           | Value                                       |
+=====================+=============================================+
|ts_matrix            | Time series matrix                          |
+---------------------+---------------------------------------------+
|wa_matrices          | List of adjecency matrices                  |
+---------------------+---------------------------------------------+
|t_lags               | Number of time lags                         |
+---------------------+---------------------------------------------+

Return Values
~~~~~~~~~~~~~~~~~~~~~~
List with lists for each spatial lag. Spatial lags lists contains the estimated spatial autocorrelation for the corresponding time lag. 
  
  List index 0 --> time lag 0 etc..
  
  
:Authors: Andreas Wolf
:Date: 2017/06/24
:Version: 1.0
