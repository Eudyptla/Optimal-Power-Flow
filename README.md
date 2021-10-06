# Introduction

Optimal power flow . 

Based on [Matpower](https://matpower.org/)

Question:

![title](image/OPF_target.jpg)

(1a) :target function

(1b) :generator active power limit.

(1c) :generator reactive power limits.

(1d) :bus voltage limits.

(1e) :swing bus angle = 0.

(1f) :active power euqation.

(1g) :reactive power euqation.

# Target Function 

1. Economics dispatch

![title](image/OPF_ED.jpg)

2. Power loss

![title](image/OPF_Ploss.jpg)

3. Reactive power preserve

![title](image/OPF_Q.jpg)

4. L-index

![title](image/OPF_L.jpg)

# Results

[9 BUS](https://matpower.org/docs/ref/matpower5.0/case9.html)

Target Function        |Cost of electricity($)|Power Loss(MW)|Reactive power(Mvar)|L-index
------------------     |:--------------------:|:------------:|:------------------:|:-------:
Original data          | 5438.3               |4.9547        |34.88               |0.1673
Economics dispatch     | 5296.7               |3.3069        |-9.611              |0.1358       
Power loss             | 6071.6               |2.3158        |-12.067             |0.1332       
Reactive power preserve| 5607.1               |2.4639        |-13.386             |0.1337         
L-index                | 6567.8               |2.3961        |-9.915              |0.1331        

[39 BUS](https://matpower.org/docs/ref/matpower5.0/case39.html)

Target Function        |Cost of electricity($)|Power Loss(MW)|Reactive power(Mvar)|L-index
------------------     |:--------------------:|:------------:|:------------------:|:--------:
Original data          | 45077                |43.641        |1274.9              |0.201
Economics dispatch     | 41857                |43.367        |1338.4              |0.19372       
Power loss             | 46815                |28.51         |1196.8              |0.1954       
Reactive power preserve| 44614                |32.578        |1095.8              |0.1935         
L-index                | 45751                |38.936        |1226.5              |0.1931   


[300 BUS](https://matpower.org/docs/ref/matpower5.0/case300.html)

Target Function        |Cost of electricity($)|Power Loss(MW)|Reactive power(Mvar)|L-index
------------------     |:--------------------:|:------------:|:------------------:|:--------:
Original data          | 72470                |408.32        |7983.7              |0.41352
Economics dispatch     | 71973                |302.82        |6971.6              |0.39503       
Power loss             | 75676                |210.6         |6195.6              |0.3974       
Reactive power preserve| 75072                |252.94        |5672.2              |0.3932         
L-index                | 73101                |364.64        |6619.2              | 0.38499   

