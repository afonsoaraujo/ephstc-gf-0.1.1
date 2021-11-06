
  real ::  Wgeq,           &! umidade do solo no ponto onde as forcas 
                 ! gravitacionais e capilares se equilibram
    Wrmax,    &! maxima umidade na cobertura vegetal
    qvsatta,        &! umidade especifica de saturacao no ar
    qvsatts,        &! umidade especifica de saturacao no solo
    qva,    &! umidade especifica do ar
    ea,   &! pressao parcial de vapor na atmosfera
    esatta,   &! pressao de saturacao no ar
    esatts,   &! pressao de saturacao no solo
    hu,             &! umidade relativa da superficie do solo
    Rg,             &! radiacao de ondas curtas que atinge o solo
    LE,             &! fluxo de calor Latente
    H,              &! fluxo de calor sensivel
    Rn,   &! fluxo de radiacao liquida
    aw,             &! calculo de qvsatta e qvsatts
    bw,             &! calculo de qvsatta e qvsatts
    Cg,   &! coeficiente termico do solo
    Ct,   &! capacidade de calor na superficie
    F1,             &!             fotossintese
    F2,             &! fatores =>  pressao hidraulica
    F3,             &!             pressao de vapor atmosferica
    F4,             &!             temperatura do ar
    f,    &! calculo de F1
    C1,             &! coeficiente do ganho de precipitacao
    C2,             &! coef. umidade contida proxima a superficie
    Pg,             &! precipitacao que alcanca o solo
    Pr,   &! precipitacao interceptada pela vegetacao
    Eg,             &! evaporacao do solo
    Etr,            &! tranpiracao da porcao "seca" das folhas
    Er,             &! evaporacao da fracao "molhada" das folhas
    Ra,   &! resistencia aerodinamica
    Rs,   &! resistencia estomatica
    Fw               ! fracao molhada do topo da vegetacao
              
  real :: G,LEg,LEr,LEtr

  real ::  Rsfc,             &! Run-off da superficie
      D12,             &! difusao entre layer 1 e 2
      D23,             &! difusao entre layer 2 e 3
      K2,             &! drenagem do layer 2 para o layer 3
      K3,             &! drenagem do layer 3 para lencol
      C3,             &! Parametros de solo
      C4,             &! Parametros de solo
      C4l,             &! Parametros de solo
      C4ref,             &! Parametros de solo
      w3_star,             &! Parametro de calculo do solo
      w23,             &! Parametro de calculo do solo
      max_d2,             &! Parametro de calculo do solo
      max_d3             ! Parametro de calculo do solo
   
  
  common /A2/ Rn,G,H,LE,LEg,LEtr,LEr,Rs,Ra
