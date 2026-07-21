This is a project I created with the PROCESSING software during my university studies.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

For those unfamiliar with PROCESSING: It is an open-source development environment used to create interactive graphics applications, data visualizations, animations, generative art, and more. It is based on the Java programming language and is designed to be easy to learn and use, even for those without extensive programming knowledge. Processing provides a wide range of functions and graphics libraries to make creating images and animations easier and faster. Additionally, third-party libraries can be integrated to extend Processing's core functionality. It was created for use by both artists and programmers and to support a wide range of applications, from digital art to scientific data processing. It is available free of charge and can run on various platforms, including Windows, Mac OS X, and Linux.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TRACK Project:

Using Processing, design an ANTHROPOMORPHIC robot operated using INVERSE KINEMATICS. The desired position (xd,yd,zd) of the gripper and its orientation must be editable via the keyboard. Specifically, for the gripper position, the xd, yd, and zd coordinates can be modified as follows: lowercase x decreases the xd coordinate, while uppercase X increases it. Similarly, the letters y and z can be used for the yd and zd coordinates.
For the desired orientation of the gripper, i.e., the Re matrix, proceed as follows. Identify the z6 axis of the gripper using the azimuth angles α and β with respect to the base system (x0,y0,z0) (refer to this figure), and define x6 and y6 followed by a rotation of an angle θ around the z6 axis. This can be demonstrated to correspond to a ZYZ-type parameterization of the desired rotation matrix Re, with angles (α, 90o-β, θ), which therefore coincides with the expression of the spherical wrist matrix R36 in which α must be substituted for θ4, 90o-β for θ5, and θ for θ6. To change the orientation of the gripper, it will therefore be sufficient to act on the variables α, β, and θ by pressing, for example, the following keys: lowercase a to decrease α and uppercase A to increase it, lowercase b to decrease β and uppercase B to increase it, lowercase t to decrease θ and uppercase T to increase it.
Write the sketch, also taking into account the following specifications:
Throughout the program's execution, the desired coordinates (xd, yd, zd) for the gripper's position and the values ​​(in degrees) of the angles (α, β, θ) that define the gripper's orientation must be displayed on the screen. The gripper's coordinates must be written using the reference triplet (x0, y0, z0) from the lesson (and not the one used by Processing).
Also display the desired Re matrix with columns of THREE DIFFERENT COLORS.
DRAW both the TRIM (x0,y0,z0) of the BASE and the TRIM (x6,y6,z6) of the CLAMP using the SAME COLORS for the three x, y, and z axes as used for the columns of Re (i.e., the x-axis should be drawn in the same color as the first column of Re, the y-axis in the same color as the second column, and the z-axis in the same color as the third column). For simplicity, the axes can be drawn without arrows.
NEGLECT the problem of COLLISIONS between the various links of the robot.
For simplicity, the CLAMP can be represented as a simple PARALLELEPIPED.
Provide the ability to switch from the HIGH ELBOW to the LOW ELBOW solution using the keyboard (for example, with the '+' and '-' keys). Instead, set the solution for the spherical wrist arbitrarily.
Include the features (already implemented in the various sketches seen in class) that allow: 1) to change the height of the view, 2) to change the value of the control law constant Kp, 3) to move the base of the robot with a mouse click.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Italian Version:

Questo é un mio progetto realizzato con il software PROCESSING durante i miei studi universitari.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Per chi non conoscesse il software PROCESSING: E'un ambiente di sviluppo open source utilizzato per creare applicazioni grafiche interattive, visualizzazioni di dati, animazioni, arte generativa e altro ancora. Si basa sul linguaggio di programmazione Java ed è progettato per essere facile da imparare e usare anche per chi non ha una conoscenza approfondita della programmazione.Processing fornisce una vasta gamma di funzioni e librerie grafiche per rendere la creazione di immagini e animazioni più facile e veloce. Inoltre, è possibile integrare librerie di terze parti per estendere le funzionalità di base di Processing.E'stato creato per essere utilizzato sia da artisti che da programmatori e per supportare una vasta gamma di applicazioni, dall'arte digitale all'elaborazione di dati scientifici. È disponibile gratuitamente e può essere eseguito su diverse piattaforme, tra cui Windows, Mac OS X e Linux.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TRACK Project:

Utilizzando Processing, disegnare un robot ANTROPOMORFO operato in CINEMATICA INVERSA. Il valore desiderato (xd,yd,zd) per la posizione della pinza e il suo orientamento devono essere modificabili da tastiera. In particolare, per quanto riguarda la posizione della pinza, le coordinate xd, yd e zd possono essere modificate nel seguente modo: con x minuscolo si diminuisce la coordinata xd, con X maiuscolo la si aumenta. Analogamente si possono usare le lettere y e z per le coordinate yd e zd.
Per quanto riguarda l'orientamento desiderato della pinza, e cioè la matrice Re, procedere come segue. Individuare l'asse z6 della pinza mediante gli angoli di azimuth α e di elevazione β rispetto al sistema di base (x0,y0,z0) (fare riferimento a questa figura), e definire x6 e y6 facendo seguire una rotazione di un angolo θ intorno all'asse z6. Questo si può dimostrare corrisponde a una parametrizzazione di tipo ZYZ della matrice di rotazione desiderata Re, con angoli (α,90o-β,θ), coincidente quindi con l'espressione della matrice R36 del polso sferico in cui occorre sostituire α al posto di θ4, 90o-β al posto di θ5 e θ al posto di θ6. Per cambiare l'orientamento della pinza sarà quindi sufficiente agire sulle variabili α, β e θ mediante la pressione per esempio dei seguenti tasti: a minuscolo per diminuire α e A maiuscolo per aumentarlo, b minuscolo per diminuire β e B maiuscolo per aumentarlo, t minuscolo per diminuire θ e T maiuscolo per aumentarlo.
Scrivere lo sketch tenendo inoltre conto delle seguenti specifiche:
Durante tutta l'esecuzione del programma deve essere riportato a schermo il valore delle coordinate desiderate (xd,yd,zd) per la posizione della pinza e quello (in gradi) degli angoli (α,β,θ) che definiscono l'orientamento della pinza. Le coordinate della pinza vanno scritte SCEGLIENDO COME TERNA DI RIFERIMENTO (x0,y0,z0) DI BASE QUELLA CONSIDERATA A LEZIONE (e non quella utilizzata da Processing).
Scrivere a schermo anche la matrice Re desiderata con le colonne di TRE COLORI DIVERSI.
DISEGNARE sia la TERNA (x0,y0,z0) della BASE sia quella (x6,y6,z6) della PINZA utilizzando per i tre assi x, y e z gli STESSI COLORI usati per le colonne di Re (cioè l'asse x va disegnato dello stesso colore della prima colonna di Re, l'asse y dello stesso colore della seconda colonna e l'asse z come la terza colonna). Per semplicità gli assi possono essere disegnati senza frecce.
TRASCURARE per semplicità il problema delle COLLISIONI tra i vari link del robot.
Per semplicità la PINZA può essere rappresentata come un semplice PARALLELEPIPEDO.
Prevedere la possibilità da tastiera (per esempio con i tasti '+' e '-') di passare dalla soluzione GOMITO ALTO a quella GOMITO BASSO. Fissare invece arbitrariamente la soluzione per il polso sferico.
Includere le funzionalità (già implementate nei vari sketch visti a lezione) che permettono: 1) di modificare l'altezza della vista, 2) di modificare il valore della costante Kp della legge di controllo, 3) di spostare la base del robot con un click di mouse.

