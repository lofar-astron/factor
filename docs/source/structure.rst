Factor Structure
================

Factor effectively sets up and runs generic pipelines that perform the actual processing. The overall structure of facet calibration as done by Factor is shown in the Figure :num:`factor-flowchart` below. The processing is divided into a number of operations, the division of which is largely determined by whether or not multiple operations may be run in parallel. In this flowchart, each operation is outlined with a black box.

.. _factor-flowchart:

.. figure:: factor_flow.pdf
   :figwidth: 90 %
   :align: center

   Factor flowchart
