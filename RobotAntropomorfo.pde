/*                                                                                                                            
                                                                                                          Progetto 2022-2023  Robot Antropomorfo che lavora in cinematica inversa Davide MARI 

Funzioni Robot Antropomorfo:
    Con le frecce Giù e Su si modifica l'altezza della vista sul robot
    Con x,y,z minuscole o maiuscole si diminuiscono o aumentano le rispettive coordinate
    Con + o - si imposta rispettivamente la configurazione gomito alto o basso
    Con le frecce Destra e Sinistra si aumenta o diminuisce la velocità Kp con cui il robot compie i propri movimenti
    
*/
                                                                                                                                    //INIZIALIZZAZIONE PARAMETRI
                                                                                                                                    
// Parametro della funzione camera() che viene modificato con le frecce SU e GIU e determina l'altezza della vista rispetto al robot
float eyeY = 0;

// Coordinate del centro del link 0 del robot che viene spostato col mouse
float xBase;
float yBase;
float zBase;

//Costante legge di controllo
float kp = 0.1;
//Step fine corsa
int stepfinecorsa=10;

// step di incremento 
float step1 = 0.1;
float step2 = 0.1;
float step3 =0.2;

// dimensioni base:
float d0_x = 250; // lungo x
float d0_y = 50; // lungo y
float d0_z = 250; // lungo z

// dimensioni link 1
float d1_x = 200; // lungo x
float d1_y = 80; // lungo y
float d1_z = 200; // lungo z

// dimensioni link 2
float d2_x = 250; // lungo x
float d2_y = 20; // lungo y
float d2_z = 20; // lungo z

// dimensioni link 3
float d3_x = 150; // lungo x
float d3_y = 20; // lungo y
float d3_z = 20; // lungo z

//dimensioni link 4
float d4_x = 100; //lungo x
float d4_y = 20; // lungo y
float d4_z = 100; // lungo z

//dimensione link 5 
float d5_x = 100; //lungo x
float d5_y = 20; //lungo y
float d5_z = 20; //lungo z

//dimensione link 6
float d6x=100; //lungo x
float d6y=20; //lungo y
float d6z=100; //lungo z

// dimensioni pinza
float d6a_x = 30; // lungo x
float d6a_y = 25; // lungo y
float d6a_z = 30; // lungo z

//dimensioni tenaglie
float d6b_x = 10; // lungo x
float d6b_y = 10; // lungo y
float d6b_z = 25; // lungo z

float posPinza = 0;
float gomito=1;

//Parametri giunto (theta1, theta2, theta3, theta4, theta5,theta6)
float[] theta = {0, 0, 0, 0, 0, 0};

//Variabili fisse(d1,d4,d6,l1,l2)
float l1 = d1_x/2;
float l2 = d2_x-d3_y/2;
float d1 = d0_y+d1_y;
float d4 = d3_x-d3_y/2+d4_x-d5_x/2;
float d6 = d5_x-d5_y/2+d6y;

//Cinematica inversa parametri necessari per trovare i vari theta successivamente
float r11;
float r21;
float r31;

float r12;
float r22;
float r32;

float r13;
float r23;
float r33;

float xd=250;
float yd=0;
float zd=250;

float x_0 = xBase; // terna di riferimento 
float y_0 = yBase; 
float z_0 = 0; 

float alpha_d=0;
float beta_d=0;
float theta_d=0;

// CILINDRO 
float bottom = 23;
float top =23; 
float h =120; // ALTEZZA POLIGONO 
int lati = 64; // NUMERI DI ANGOLI DEL POLIGONO 

void setup()
{
  size(1000, 800, P3D);
  stroke(50);
  strokeWeight(2);
  xBase = width/2;
  yBase = height/2;
  
}


void draw()
{
  background(0);
  directionalLight(128, 128, 128, -1, 0, 1);
  lights();
  // Permette di ruotare la vista:
  camera((width/2.0), height/2 - eyeY, (height/2.0) / tan(PI*60.0 / 360.0), width/2.0, height/2.0, 0, 0, 1, 0);

  if(mousePressed) {
    xBase=mouseX;
    yBase=mouseY;
  }
  if(keyPressed)
  {
    //movimeto camera
    if (keyCode==DOWN){
    eyeY -= 5;
    }
    if(keyCode == UP){
    eyeY += 5;
    }
    if(key == '+'){
      gomito = 1;  //gomito alto
    }
    if(key == '-'){
      gomito = -1; //gomito basso
    }
    if (key == 'x'){
      xd -= stepfinecorsa*kp;
    }
    if (key == 'X'){
      xd += stepfinecorsa*kp;
    }
    if (key == 'y'){
      yd -= stepfinecorsa*kp;
    }
    if (key == 'Y'){
      yd += stepfinecorsa*kp;
    }
    if (key == 'z'){
      zd -= stepfinecorsa*kp;
    }
    if (key == 'Z'){
      zd += stepfinecorsa*kp;
    }
    if (key == 'a'){
      alpha_d -= step1*kp;
    }
    if (key == 'A'){
      alpha_d += step1*kp;
    }
    if(key == 'b'){
      beta_d -= step2*kp;
    }
    if(key == 'B'){
      beta_d+= step2*kp;
    }
    if(key == 't'){
      theta_d -= step3*kp;
    }
    if(key == 'T'){
      theta_d += step3*kp;
    }
    if(key == 'K'){
      kp += 0.01;
      if(kp>0.5){
        kp = 0.5;
      }
    }
    if(key == 'k'){
      kp -= 0.01;
      if(kp<0.1){
        kp=0.1;
      }
    }
  }
  Testi();
                                                                                                                                          //Inizializzo Matrice Re
//Inizializzo delle variabili che mi servono per creare la matrice Re richiesta nel testo(lo faccio per non ripetere i vari sin e cos al momento della creazione di Re)                                                                                                                                //
      r11= cos(alpha_d)*sin(beta_d)*cos(theta_d)-sin(alpha_d)*sin(theta_d); 
      r21= sin(alpha_d)*sin(beta_d)*cos(theta_d)+cos(alpha_d)*sin(theta_d); 
      r31= -cos(beta_d)*cos(theta_d);
      
      r12= -cos(alpha_d)*sin(beta_d)*sin(theta_d)-sin(alpha_d)*cos(theta_d); 
      r22= -sin(alpha_d)*sin(beta_d)*sin(theta_d)+cos(alpha_d)*cos(theta_d);
      r32= cos(beta_d)*sin(theta_d);
      
      r13= cos(alpha_d)*cos(beta_d);
      r23= sin(alpha_d)*cos(beta_d);
      r33= sin(beta_d);
      
float [][] Re = {{r11,r12,r13},
                {r21,r22,r23},
                {r31,r32,r33}};  
   
//Funzione che gestisce la cinematica inversa
  cinematicainversa(Re);
  
// scrivo la matrice di rotazione utilizzando la funzione creata
  
scriviMatrice("Re :", Re, 20,350);

//Funzione per gestione Robot Antropomorfo
    ANTROPOMORFO();
 
}
                                                                                                    //FUNZIONE PER LA CREAZIONE DEL ROBOT ANTROPOMORFO
void ANTROPOMORFO() {
  
 fill(#00FFFF);
 strokeWeight(2);
 translate(xBase,yBase,zBase);  
 rotateX(PI/2);
 riferimento_iniziale(); // disegno sistema di riferimento iniziale 
 strokeWeight(2);
 stroke(255);
 box(d0_x, d0_z, d0_y);
 
  //link1
  fill(#F07400);
  translate(0,0,d0_y/2);
  rotateZ(-theta[0]);
  translate(0,0,d1_y/2);
  box(d1_x,d1_z,d1_y);
  
 // Link 2 (si muove con THETA 2 = theta[1]) 
  translate(d1_x/2,0, d1_y/2);
  rotateX(-PI/2); //giunto2:z uscente, y giù e x a destra
  strokeWeight(2);
  stroke(255);
  rotateZ(-theta[1]);
  cilindro(); // disegno cilindri su i giunti rotoidali 
  translate(d2_x/2,0,0);
  box(d2_x,d2_y,d2_z); // crea il poligono

 // Link 3 (si muove con THETA 3 = theta[2])
  translate((d2_x-d3_y)/2, 0, 0);
  rotateZ(-theta[2]); //giunto3=giunto2
  rotateX(-PI/2);//giunto4:z giù, y entrante e x a destra
  strokeWeight(2);
  stroke(255);
  translate(0, 0, (d3_x-d3_y)/2);
  box(d3_z,d3_y,d3_x); 

 // Link 4 (si muove con THETA 4 = theta[3])
 translate(0, 0, d3_x/2);
  rotateZ(-theta[3]);
  translate(0, 0, d4_x/2);
  box(d4_z,d4_y,d4_x);
  
  //link5
  translate(0, 0, d4_x/2-d5_z/2);
  //Per quanto concerne giunto 5 z uscente,y verso giù e x verso destra
  rotateX(PI/2);
  //z uscente, y verso su,x verso sinistra
  rotateZ(PI);
  strokeWeight(2);
  stroke(255);
  rotateZ(-theta[4]);

  //giunto 5
  translate(0, -(d5_x/2-d5_z/2), 0);
  box(d5_y, d5_x, d5_z);
  
  //link6
  translate(0,-d5_x/2,0);
  //Per quanto concerne giunto 6 z verso giù,y uscente e x verso sinistra
  rotateX(PI/2);
  strokeWeight(2);
  stroke(255);
  rotateZ(-theta[5]); 
  translate(0, 0, d6y/2);
  box(d6x,d6z,d6y);
  
  translate(0, 0, d6y/2);
  rotateZ(PI/2);
  riferimento_finale(); // disegno sistema di riferimento finale
  strokeWeight(2);
  stroke(255);
  Pinza(); // disegno pinza 
}  

void Pinza(){

     // Pinza FORCHE PINZA
  pushMatrix();
  translate(-d6b_x/2-posPinza, d6y/2+d6b_y/2, 0);
  box(d6b_x,d6b_y,d6b_z);
  popMatrix();
  translate(d6b_x/2+posPinza, d6y/2+d6b_y/2,0);
  box(d6b_x,d6b_y,d6b_z);
  
}

void riferimento_iniziale()
{
  
 //Sistema di riferimento della base, link 0  
   strokeWeight(4);   
   stroke(255,0,0);    //asse x rosso
   line(0, 0, 0,200, 0, 0);
   text("X_0", 210, 0, 0);
   stroke(0,255,0);  //asse y verde
   line(0,0,0, 0, -200, 0);
   text("Y_0", 0, -210, 0);
   stroke(0,0,255); //asse z blu
   line(0,0,0,0, 0, 120);
   text("Z_0", 0, 0,130);
}

void riferimento_finale()
{
  //Sistema di riferimento link6
   strokeWeight(4);
   stroke(255,0,0);    //asse x rosso
   line(0, 0, 0,0, 100, 0);
   text("X_6", 0, 110, 0);
   stroke(0,255,0);  //asse y verde
   line(0,0,0, 200, 0, 0);
   text("Y_6", 210, 0, 0);
   stroke(0,0,255); //asse z blu
   line(0,0,0,0, 0,120);
   text("Z_6", 0, 0,130);
}

                                                                    // FUNZIONE PER LA CINEMATICA INVERSA(Uso formule studiate inerenti alla cinematica inversa del robot antropomorfo)

void cinematicainversa(float [][]Re)
{
//Variabili che mi servono per calcolo dei theta
  float Valore1;
  float Valore2;
  float Valx=xd-(d6z*Re[0][2]);
  float Valy=yd-(d6z*Re[1][2]);
  float Valz=zd-(d6z*Re[2][2]);
  float argAcos;
  
  //calcolo theta1
  theta[0]=atan2(Valy,Valx);
  
  //calcolo theta3
 Valore1=Valx*cos(theta[0])+Valy*sin(theta[0])-l1;
 Valore2=d1-Valz;
  
  argAcos=(Valore1*Valore1+Valore2*Valore2-d4*d4-l2*l2)/(2*l2*d4);
  if (argAcos < -1 || argAcos > 1)
  {
    textSize(25);
    fill(256, 0, 0);
    text("POSIZIONE FUORI DALLO SPAZIO DI LAVORO",   800, 800);
  } 
   else {
      if(gomito == 1){
        theta[2]=asin(argAcos);
      }
      else {
        theta[2]=PI-asin(argAcos);
      }

  }
  //calcolo theta2
  
  theta[1] = atan2(d4*cos(theta[2])*Valore1-(d4*sin(theta[2])+l2)*Valore2 , (d4*sin(theta[2])+l2)*Valore1+d4*cos(theta[2])*Valore2);
  
  //Calcolo matrici R36 e R03 che utilizzero per calcolo degli ultimi theta
  
  float R03[][] = {{cos(theta[0])*cos(theta[1]+theta[2]), sin(theta[0]), cos(theta[0])*sin(theta[1]+theta[2])},
                   {sin(theta[0])*cos(theta[1]+theta[2]), -cos(theta[0]), sin(theta[0])*sin(theta[1]+theta[2])},
                   {sin(theta[1]+theta[2]), 0, -cos(theta[1]+theta[2])}};
 //UTILIZZO FUNZIONE AUSILIARA PER TROVARE LA TRASPOSTA          
  float R03_t [][] = trasposta(R03);
  
  float R36 [][] = mProd(R03_t,Re);
  
 //Utilizzo funzione ausiliaria che mi peremette di trascrivere a schermo una matrice qualsiasi con i valori passati
  scriviMatrice("R36: ", R36, 20, 550);
  scriviMatrice("R03: ", R03, 20, 800);
  
  
  float theta5=acos(R36[2][2]);
  //verifico che theta5 sia diverso da 0 e +-180 poichè dovrò usare formule diverse per trovare i theta restanti
  if(theta5!=0 && theta5!=PI && theta5!=-PI) { 
    if(sin(theta5)>0) { //theta5 appartiene all'intervallo (0,180)
    
      //calcolo theta5
      theta[4] = atan2(sqrt(R36[0][2]*R36[0][2]+R36[1][2]*R36[1][2]),R36[2][2]);
      
      //calcolo theta4
      theta[3] = atan2(R36[1][2],R36[0][2]);
      
      //calcolo theta6
      theta[5] = atan2(R36[2][1], -R36[2][0]);
    }
    else if(sin(theta5)<0) {
      
      //calcolo theta5
      theta[4] = atan2(-sqrt(R36[0][2]*R36[0][2]+R36[1][2]*R36[1][2]),R36[2][2]);
      
      //calcolo theta4
      theta[3] = atan2(-R36[1][2],-R36[0][2]);
      
      //calcolo theta6
      theta[5] = atan2(-R36[2][1], +R36[2][0]);
    }
    else{
      if (theta5 == 0){
      theta[5]=PI;
      theta[3]=atan2(R36[0][1],R36[0][0])-theta[5];
    }else if(theta5==PI || theta5==-PI){
      theta[5]=PI;
      theta[3]=atan2(-R36[0][1],-R36[0][0])+theta[5];
    }
    }
  }
  
} 
void cilindro()
{
  pushMatrix();
  noStroke();
  rotateX(PI/2);
 
  float angoli;
  float x[] = new float[lati+1];
  float z[] = new float[lati+1];
  
  float x2[] = new float[lati+1];
  float z2[] = new float[lati+1];
 
  //get the x and z position on a circle for all the sides
  for(int i=0; i < x.length; i++){
    angoli = 2*PI / (lati) * i;
    x[i] = sin(angoli) * bottom;
    z[i] = cos(angoli) * bottom;
  }
  
  for(int i=0; i < x.length; i++){
    angoli = TWO_PI / (lati) * i;
    x2[i] = sin(angoli) * top;
    z2[i] = cos(angoli) * top;
  }
 
  //draw the bottom of the cylinder
  beginShape(TRIANGLE_FAN);
 
  vertex(0,   -h/2,    0);
 
  for(int i=0; i < x.length; i++){
    vertex(x[i], -h/2, z[i]);
  }
 
  endShape();
 
  //draw the center of the cylinder
  beginShape(QUAD_STRIP); 
 
  for(int i=0; i < x.length; i++){
    vertex(x[i], -h/2, z[i]);
    vertex(x2[i], h/2, z2[i]);
  }
 
  endShape();
 
  //draw the top of the cylinder
  beginShape(TRIANGLE_FAN); 
 
  vertex(0,   h/2,    0);
 
  for(int i=0; i < x.length; i++){
    vertex(x2[i], h/2, z2[i]);
  }
 
  endShape();
  stroke(255);
  strokeWeight( 1);
  popMatrix();
}

void Testi()
{
  //Inserimento testi 
  textSize(25);
  fill(255,0,0);
  
  //KP di controllo;
  text("kp = ",10, 30);
  text(kp,60,30);
  
  // terna di riferimento 

  text("x_0 = ",20,100); 
  text(x_0,80,100);
  
  text("y_0 = ",20,125); 
  text(y_0,80,125);
 
  text("z_0 = ",20,150); 
  text(z_0,80,150);
  
  fill(0,255,0);  
  text("coordinata y vista = ",600,30); 
  text(eyeY,900,30);
  text("percentuale apertura pinza = ",600,80); 
  text(round(100*posPinza/(d6a_x/2-d6b_x)),980,80);
  text("%",1025,80);
  
  text("Alpha_d :", 1250, 30);
  text(alpha_d*180/PI, 1350, 30);
  text("Beta_d  :", 1250, 50);
  text(beta_d*180/PI, 1350, 50);
  text("Theta_d :", 1250, 70);
  text(theta_d*180/PI, 1350, 70);

  textSize(25);
  fill(255,0,0);
  
  text("xd =  ",320,20); 
  text(xd,430,20);
  
  text("yd = ",320,70); 
  text(yd,430,70);
  
  text("z_d = ",320,120); 
  text(zd,430,120);
  
}

                              ///////////////////////////////////////////////////////////////////////////////////FUNZIONI AUSILIARIE /////////////////////////////////////////////////////////////////
float[][] mProd(float[][] A,float[][] B) // Calcola prodotto di due matrici A e B
{
  int nA = A.length;
  int nAB = A[0].length;
  int nB = B[0].length;
  
  float[][] C = new float[nA][nB]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nB; j++) 
    {  
      for (int k=0; k < nAB; k++) 
      {
        C[i][j] += A[i][k] * B[k][j];
      }
    }
  }
  return C;
}

float[][] mSum(float[][] A,float[][] B) // Calcola la somma di due matrici A e B
{
  int nA = A.length;
  int nB = A[0].length;
  
  float[][] C = new float[nA][nB]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nB; j++) 
    {  
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return C;
}

float[][] trasposta(float[][] A) // Calcola la trasposta di una matrice A
{
  int nR = A.length;
  int nC = A[0].length; 
  
  float[][] C = new float[nC][nR]; 

  for (int i=0; i < nC; i++) 
  {
    for (int j=0; j < nR; j++) 
    {  
      C[i][j] = A[j][i];
    }
  }
  return C;
}


float[][] minore(float[][] A, int i, int j) // Determina il minore (i,j) di una matrice A
{
  int nA = A.length;
  float[][] C = new float[nA-1][nA-1];
  
  for (int iM = 0; iM < i; iM++)
  {
    for (int jM = 0; jM < j; jM++)
    {
      C[iM][jM] = A[iM][jM];
    } 
    for (int jM = j; jM < nA-1; jM++)
    {
      C[iM][jM] = A[iM][jM+1];
    } 
  }
  for (int iM = i; iM < nA-1; iM++)
  {
    for (int jM = 0; jM < j; jM++)
    {
      C[iM][jM] = A[iM+1][jM];
    } 
    for (int jM = j; jM < nA-1; jM++)
    {
      C[iM][jM] = A[iM+1][jM+1];
    } 
  }
  return C;
}


float det(float[][] A) // Calcola il determinante di A
{
  int nA = A.length;
  float determinante = 0;
  
  if (nA == 1)
  {
    determinante = A[0][0];
  }
  else
  {
    for (int j=0; j < nA; j++) 
    {
      determinante = determinante + A[0][j]*pow(-1,j)*det(minore(A,0,j));
    }
  }
  return determinante;
}


float[][] invMat(float[][] A) // Calcola l'inversa di una matrice A
{
  int nA = A.length;
  float[][] C = new float[nA][nA];
  float detA = det(A);

  if (nA == 1)
  {
    C[0][0] = 1/detA;
  }
  else
  {
    for (int i=0; i < nA; i++) 
    {
      for (int j=0; j < nA; j++) 
      {
        C[j][i] = pow(-1,i+j)*det(minore(A,i,j))/detA;
      }
    }
  }
  return C;
}

float[][] idMat(int nA, float sigma) // Assegna una matrice identità di ordine nA moltiplicata per una costante sigma
{
  float[][] I = new float[nA][nA]; 

  for (int i=0; i < nA; i++) 
  {
    for (int j=0; j < nA; j++) 
    {  
      I[i][j] = 0;
    }
    I[i][i] = sigma;
  }
  return I;
}

void scriviMatrice(String s, float[][] M, int x, int y) // Scrive una matrice a partire dal punto (x,y)
{
  textSize(20);
  fill(255);
  text(s,x,y);
  fill(255,0,0);
  text(M[0][0],x,y+30); text(M[1][0],x,y+60); text(M[2][0],x,y+90);              
  fill(0,255,0);
  text(M[0][1],x+90,y+30); text(M[1][1],x+90,y+60); text(M[2][1],x+90,y+90);                               
  fill(0,0,255);
  text(M[0][2],x+180,y+30); text(M[1][2],x+180,y+60); text(M[2][2],x+180,y+90);
}
  
