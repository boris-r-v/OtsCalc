package ots_calc;
/** Класс для расчёта мгновенной схемы ОТС
 на постоянном токе НЕ однородный рельс
 */
public class OTS_DC_HeterogRail_call {
  /* В комментариях в данном классе используются следующие обозначения:
  Р - рельсы данного пути, рассматривается как длинная линия
  ЭПС - ЭПС по данному пути
  ФОТ - фидера обратного тока (отсасывающие фидера) подключенные к рельсам данного пути
  ЗАЗ - заземлители к рельсам данного пути (имеет смысл учитывать заземлители с малым сопротивлением заземления менее 10 Ом)
  МПС -  междупутные соединители, условное название – фактически это сопротивление, соединяющее две любые точки каких любо путей (длинных линий)
  МКР – метод конечных разностей, сеточный метод дискретизации дифференциальных уравнений (в данном случае уравнение длинной электрической линии) по пространству (здесь вдоль пути)


  Единицы измерения в классе:
  координаты и длины - км;
  напряжение  - В;
  ток - А;
  сопротивление - Ом;
  проводимость - См;
  погонное сопротивление - Ом/км;
  переходное сопротивление Ом*км;
  */


    //-----------------------------------------------------Головной класс переменные и инициализация-------------------------------------------------------------------------

    // переменные дочерних классов
    public Computing_settings computing_settings;
    public Errors_and_messages err_and_mes;
    public Mesh[] meshes; // массив сеток
    public Track[] tracks; // массивы путей

    //переменные стандартных типов данных
    public double mps[][] = {{}}; // двумерный массив в каждой строке 5 элементов: 0 - номер пути начальной точки подключения, 1 - номер пути конечной точки подключения, 2 - координата начальной точки подключения км, 3 -координата конечной точки подключения км, 4 - сопротивление Ом.
    private double I_mps[]; // массив токов междупутных соединителей МПС и токов отходящих ветвей в местах соединения с главными путями
    private int Ntrc; // количество главных путей и ответвлений
    private int Nmsh; // количество сеток
    private int num_mesh_mps[][], // двумерный массивы номеров сеток для МПС (в каждой строке начальная и конечная точка)
            num_track_mps[][], // двумерный номеров путей для МПС (в каждой строке начальная и конечная точка)
            index_mps[][]; //двумерный массивы  индексов узлов по сетке МПС (в каждой строке начальная и конечная точка)
    private double a_x_find[][]; // матрица коэффициентов влияния тока во всех МПС на напряжения во всех МПС Ом. По главной диагонали сами на себя.
    private double U_const[]; //  массивы разности напряжений от заданных токов в МПС от начальнйо до конечной точки подключения

    // инициализатор основного класса по умолчанию
    OTS_DC_HeterogRail_call(int Ntrc, int Nmsh) { // количество путей и сеток
        this.Ntrc = Ntrc;
        this.Nmsh = Nmsh;
        this.computing_settings = new Computing_settings();
        this.err_and_mes = new Errors_and_messages();
        this.meshes = new Mesh[Nmsh];
        this.tracks = new Track[Ntrc];
        for (int i = 0; i < Ntrc; ++i) { // заполняем массив путей
            this.tracks[i] = new Track();
        }
        for (int i = 0; i < Nmsh; ++i) { // заполняем массив сеток
            this.meshes[i] = new Mesh();
        }
    }


    //-----------------------------------------------------Подклассы, их методы, свойства и инициализаторы-------------------------------------------------------------------------

    /*совмещенный класс для пути
     */
    public class Track {
        public double r[][], //функция погонного сопротивления рельсов вдоль пути. В каждой строке два значения: 1- координата км, 2 -величина погонного сопротивление рельсов Ом/км до данной координаты.
                rp[][], // функция переходного сопротивления рельсы-земля вдоль пути. В каждой строке два значения: 1- координата км, 2 -величина переходного сопротивления Ом*км до данной координаты.
                fot[][], // таблица подключенных к данному пути ФОТ. В каждой строке два значения: 1- координата точки подключения км, 2 -величина тока ФОТ А. Положительный ток принят при работе ТП в режиме выпрямителя (когда ток из рельсов затекает в ФОТ)
                eps[][], // таблица ЭПС находящихся на данном пути. В каждой строке два значения: 1- координата расположения км, 2 -величина тока ЭПС А. Положительный ток принят при работе ЭПС в режиме тяги или хх (когда ток стекает из колесных пар в рельсы)
                zaz[][], // таблица ЗАЗ находящихся на данном пути. В каждой строке два значения: 1- координата точки подключения км, 2 -величина сопротивления заземления Ом.
                R_tch[][], // таблица сосредоточенных точечных сопротивлений в рельсах на данном пути. В каждой строке два значения: 1- координата расположения км, 2 -величина сосредоточенного сопротивления Ом. (к примеру: разрыв рельсовой сети, дефектный неизолированный стык с высоким сопротивлением, дефект в средней точке ДТ с повышенным сопротивлением).
                Rv0, // волновое сопротивление слева Ом
                Rvn; // волновое сопротивление справа Ом
        public int num_mesh; // номер сетки для данного пути
        private double m3db[][], // матрица 3диагональная ленточная См
                vector_b[], // вектор правой части А
                u[], // искомый вектор напряжения в узлах сетки  В
                i[], // вектор тока рельсов в узлах сетки А
                i_grnd[]; // вектор тока земле в узлах сетки от данного пути А

        Track() {
            fot = new double[][]{};
            eps = new double[][]{};
            zaz = new double[][]{};
            R_tch = new double[][]{};
            Rv0 = -1; // присвоим волновое сопротивление в начале и в конце значение -1, что автоматически означает расчёт по величине r и rp
            Rvn = -1;
        }
    }


    /** Класс расчётная сетка
     дискретизация вдоль пути рельсов как электрически длинной линии
     линейная однородная сетка как основа метода конечных разностей
     */
    public class Mesh {
        public double X_beg, X_end; //начальная и конечная координата участка расчёта км
        public double dX; //шаг сетки по длине в км
        private double X[]; // массив координат узлов сетки
        private int mesh_N; // количество узлов сетки

        Mesh() {dX = 0.1;
            X_beg = 0;
            X_end = 20; } // по умолчанию параметры сетки такие

        /*метод создает сетку
        за основу берет данные из полей класса: X_beg,X_end, dX
        */
        private void create_mesh() {
            double X0 = dX*(int) Math.round(X_beg/dX);
            double Xn = dX*(int) Math.round(X_end/dX);
            this.mesh_N = (int) ((Xn-X0)/dX)+1;
            X = new double[this.mesh_N];
            X[0] = X0;
            for (int i = 1; i < this.mesh_N; ++i) {
                X[i] = X[i-1]+dX;
            }
        }

        /*метод распределяет функцию заданную таблицей table[][] на сетку
        table[][] это таблица в каждой строке которой содержится два элемента координата км и значение функции
        координата в каждой следующей строке должна возрастать
         table[][] должен содержать как минимум одну строку, последняя строка считается до конца сетки не зависимо от координаты в ней
        */
        private double[] distribute_function_over_mesh(double table[][]) {
            int N = this.X.length;
            double out[] = new double[N];
            if (table.length == 1) { // если в таблице только одна строка, то по всей сетке одно значение
                for (int i = 0; i < N; ++i) {
                    out[i] = table[0][1];
                }
                return out;
            }
            //если как минимум две строки
            int N_index = table.length-1, k = 0,
                    indexes[] = new int[N_index]; // создадим массив индексов по сетке для каждой строки за исключением последней
            for (int i = 0; i < N_index; ++i) {
                indexes[i] = find_near_index_over_mesh(table[i][0]);
                if (indexes[i] >= 0) {
                    while (k <= indexes[i]) {
                        out[k] = table[i][1];
                        k++;
                    }
                }
            }
            while (k < N) {
                out[k] = table[N_index][1];
                k++;
            }
            return out;
        }

        // находит ближайший индекс элемента сетки по заданной координате  X
        private int find_near_index_over_mesh(double X) {
            if (X < this.X[0]) { // если выходит за левую границу возвращает -1
                return -1;
            }
            if (X > this.X[this.X.length-1]) { // если выходит за правую границу возвращает -2
                return -2;
            }
            return (int) Math.round((X-this.X[0])/this.dX);
        }

        /* находит два ближайших индекса элемента сетки по заданной координате  X
         возвращает int[2] - номер ближайшего элемента слева и справа
         */
        private int[] find_2near_index_over_mesh(double X) {
            if (X < this.X[0]) { // если выходит за левую границу
                return new int[]{-1, -1};
            }
            if (X > this.X[this.X.length-1]) { // если выходит за правую границу
                return new int[]{-2, -2};
            }
            return new int[]{(int) Math.floor((X-this.X[0])/this.dX), (int) Math.ceil((X-this.X[0])/this.dX)};
        }

        /*создание 3Диганальной матрицы в ленточном виде
        это матрица  состоит из трех строк: нижняя диагональ, главная диагональ и верхняя диагональ
        на вход принимает экземпляр вложенного класса Track
      */
        private void create_3diag_matrix_band(Track track) {
            int N = this.mesh_N, index;
            double r[] = distribute_function_over_mesh(track.r), rp[] = distribute_function_over_mesh(track.rp), // Распределение по сетке функций сопротивления рельса Ом/км и переходное сопротивление рельс-земля Ом*км
                    diag_dw[] = new double[N], diag[] = new double[N], diag_up[] = new double[N], // три вектора: нижняя (под главной) диагональ, главная диагональ, верхняя (над главной диагональю). Эти вектора будут иметь размерность проводимости См
                    dX = this.dX, // шаг по сетке [км]
                    r0, // волновое сопротивление Rv0 [Ом] и сопротивление r0 [Ом] элемента рельса в начале
                    rn_mn, // тоже самое в конце
                    r_i_mn, r_i, rp_i; // сопротивление i-ого элемента рельса и его сопротивление заземления все в Ом

            //формируем массив точечных сосредоточенных сопротивлений в рельсах
            double R_tch[] = new double[N-1];
            for (int i = 0; i < track.R_tch.length; ++i) {
                index = this.find_2near_index_over_mesh(track.R_tch[i][0])[0];
                R_tch[index] += track.R_tch[i][1];
            }
            // далее все распределенные параметры дискретезуются по сетке в сосредоточенные для каждого элемента
            //из Ом/км -> Ом (вдоль пути); Ом*км -> Ом (на землю)
            if (track.Rv0 == -1) { // если волновое сопротивление =-1, то определим автоматически
                track.Rv0 = Math.sqrt(r[0]*rp[0]); // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
            }
            if (track.Rvn == -1) {
                track.Rvn = Math.sqrt(r[N-1]*rp[N-1]); // sqrt(Ом/км*Ом*км)=sqrt(Ом*Ом)=Ом
            }
            r0 = 0.5*(r[0]+r[1])*dX+R_tch[0]; // Ом/км*км=Ом
            rn_mn = 0.5*(r[N-2]+r[N-1])*dX+R_tch[N-2];

            //заполняем первый столбец трех диагоналей
            //все элементы трехдиагональной матрицы имеют размерность См=1/Ом
            diag_up[0] = 0;
            diag[0] = 1/track.Rv0+1/r0+0.5*dX/rp[0]; //  1/Ом+1/Ом+км/(Ом*км)=См
            diag_dw[0] = -1/r0; //  1/Ом
            //заполняем последний столбец трех диагоналей
            diag_up[N-1] = -1/rn_mn; //  1/Ом=См
            diag[N-1] = 1/track.Rvn+1/rn_mn+0.5*dX/rp[N-1]; // 1/Ом+1/Ом+км/(Ом*км)=См
            diag_dw[N-1] = 0;
            //заполянем остальные столбцы
            for (int i = 1; i < N-1; ++i) {
                r_i_mn = 0.5*(r[i]+r[i-1])*dX+R_tch[i-1]; // Ом/км*км=Ом
                r_i = 0.5*(r[i]+r[i+1])*dX+R_tch[i]; // Ом/км*км=Ом
                rp_i = rp[i]/dX; // Ом*км/км=Ом
                diag[i] = 1/r_i_mn+1/r_i+1/rp_i; //заполняем главную диагональ 1/Ом=См
                diag_dw[i-1] = -1/r_i_mn; //заполняем нижнюю диагональ
                diag_up[i+1] = -1/r_i; //заполняем верхнюю диагональ
            }
            diag_up[1] = -1/r0; // дозаполняем нижнюю и верхнюю диагональ строка 2 и предпоследняя
            diag_dw[N-2] = -1/rn_mn;

            // добавляем проводимости заземлителей zaz в главную диагональ
            for (int i = 0; i < track.zaz.length; ++i) {
                index = this.find_near_index_over_mesh(track.zaz[i][0]);
                diag[index] += 1/track.zaz[i][1];
            }
            track.m3db = new double[][]{diag_dw, diag, diag_up};
        }

        /*геттер возвращает массив координат узлов сетки
         */
        public double[] get_X() {
            return this.X;
        }
    }

    /** класс вычислительных настроек для электрического расчёта ОТС
     */
    class Computing_settings {
        public double convergence_U; // допустимая невязка величины напряжения усреднённая по всем точкам граничного условия, в которых происходит итерационный поиск величины тока втекающего в эти точки.
        public int max_number_iterat; // максимальное число итераций при расчёте величины тока, втекающего в граничные точки, по методу Ньютона в цикле.
        public double initial_damping_factor; // начальное коэффициент демпфирования в методе Ньютона при расчёте величины тока, втекающего в граничные точки.
        private double current_state_solver[]; // текущее состояние решателя содержит массив из трёх чисел полученных в конце последнего расчёта: 1 - количество итераций, 2- средняя невязка по напряжению, 3 - коэффициент демпфирования

        Computing_settings() // инициализатор по умолчанию
        {convergence_U = 0.01; initial_damping_factor = 0.7; max_number_iterat = 1000; }

        //возвращает параметры решаттеля в виде массива
        public double[] get_current_state_solver() {
            return current_state_solver;
        }
    }

    /** класс для ошибок при вычислениях, вводе данных, и сообщений о них
     */
    class Errors_and_messages {
        private boolean data_error, solver_error, calc_completed; //data_error - ошибка в данных (к примеру, неправильная длина массива), solver_error - шибка решателя (к примеру, не достигнута сходимость, исчерпано число итераций), calc_completed  - признак что расчёт выполнен хотя бы раз
        private String messeg_data_error, messeg_solver_error; // текстовое сообщение об этих ошибках
        Errors_and_messages() // инициализатор по умолчанию - ошибок нет, расчёт не выполнен
        {data_error = false; solver_error = false; calc_completed = false;
            messeg_data_error = ""; messeg_solver_error = "Расчёт не выполнен"; }
        // геттеры для всех свойств
        public boolean get_data_error() {return data_error; }
        public boolean get_solver_error() {return solver_error; }
        public String get_messeg_data_error() {return messeg_data_error; }
        public String get_messeg_solver_error() {return messeg_solver_error; }

        //обнуление для ошибки данных
        private void reset_data_error() {
            data_error = false; messeg_data_error = "";
        }
        //обнуление для ошибки решателя
        private void reset_solver_error() {
            solver_error = false; messeg_solver_error = "";
        }
    }


    //-----------------------------------------------------Внутренние процедуры базового класса приватные и публичные (за исключением сеттеров и геттеров)-------------------------------------------------------------------------

    //проверка исходных данных
    private void verify_data(boolean verify_I, boolean verify_Xfot_mps, boolean verify_Xeps) { // verify_I - условие проверки токов; verify_Xfot_mps - проверка коорд ФОТ, МПС, ЗАЗ; verify_Xeps - проверка коорд ЭПС
        err_and_mes.reset_data_error(); // очистка data_error

        // проверка по заданию ФОТ и ЭПС и токам
        int Neps_fot_all = 0;
        double Isum_fot = 0, Isum_eps = 0, deltaI;
        if (verify_I) {
            // проверяем чтобы сумма токов ЭПС равнялась сумме токов ФОТ не превышая погрешности
            // и что вообще ЭПС  и ФОТ заданы
            for (int i = 0; i < this.Ntrc; ++i) {
                Neps_fot_all += tracks[i].fot.length+tracks[i].eps.length;
                for (int j = 0; j < tracks[i].fot.length; ++j) {
                    Isum_fot += tracks[i].fot[j][1];
                }
                for (int j = 0; j < tracks[i].eps.length; ++j) {
                    Isum_eps += tracks[i].eps[j][1];
                }
            }
            if (Neps_fot_all == 0) { // если массивы ФОТ и/или ЭПС не заданы
                err_and_mes.data_error = true;
                err_and_mes.messeg_data_error += "Массивы ФОТ и/или ЭПС не заданы. ";
            }
            deltaI = Math.abs(Isum_fot-Isum_eps);
            if (deltaI/(Math.max(Math.abs(Isum_fot), Math.abs(Isum_eps))+1e-6) > 0.05) {
                err_and_mes.messeg_data_error += " Предупреждение: сумма токов ФОТ и ЭПС расходится более чем на 5 %";
            }
        }

        // проверки по сетке и координатам
        double all_points[][][];
        int num_mesh1, num_mesh2,
                num_track1, num_track2;
        if (verify_Xfot_mps) {
            // проверяем чтобы  границы сетки заданы корректно
            for (int i = 0; i < this.Nmsh; ++i) {
                if ((this.meshes[i].X_end-this.meshes[i].X_beg) <= 2*this.meshes[i].dX) { // если ошибка в границах сетки и менее двух узлов
                    err_and_mes.data_error = true;
                    err_and_mes.messeg_data_error += "Сетка номер "+i+": границы сетки заданы не корректно, либо получается менее трёх узлов";
                }
            }
            //проверяем чтобы координаты точек ФОТ, ЗАЗ, R_tch укладывались в границы сетки данного пути
            for (int i = 0; i < this.Ntrc; ++i) {
                all_points = new double[][][]{tracks[i].fot, tracks[i].zaz, tracks[i].R_tch};
                num_mesh1 = tracks[i].num_mesh;
                for (int j = 0; j < all_points.length; ++j) {
                    for (int k = 0; k < all_points[j].length; ++k) {
                        if ((all_points[j][k][0] > this.meshes[num_mesh1].X_end) || (all_points[j][k][0] < this.meshes[num_mesh1].X_beg)) {
                            err_and_mes.data_error = true;
                            err_and_mes.messeg_data_error += "Путь номер "+(i+1)+": координаты точек ФОТ, ЗАЗ или сосред сопротивлен выходят за границы сетки";
                        }
                    }
                }
            }
            // проверка чтобы точки подключения МПС к путям в пределах сетки
            for (int i = 0; i < this.mps.length; ++i) {
                num_track1 = (int) this.mps[i][0];
                num_track2 = (int) this.mps[i][1];
                num_mesh1 = tracks[num_track1].num_mesh;
                num_mesh2 = tracks[num_track2].num_mesh;
                if ((this.mps[i][2] > this.meshes[num_mesh1].X_end) ||
                        (this.mps[i][2] < this.meshes[num_mesh1].X_beg)) {
                    err_and_mes.data_error = true;
                    err_and_mes.messeg_data_error += "МПС "+(i)+": координата начальной точки подключения к пути "+num_track1+"выходbт за границы сетки";
                };
                if ((this.mps[i][3] > this.meshes[num_mesh2].X_end) ||
                        (this.mps[i][3] < this.meshes[num_mesh2].X_beg)) {
                    err_and_mes.data_error = true;
                    err_and_mes.messeg_data_error += "МПС "+(i)+": координата конечной точки подключения к пути "+num_track2+"выходbт за границы сетки";
                };
            }
        }

        if (verify_Xfot_mps) {
            //проверяем чтобы координаты точек ЭПС укладывались в границы сетки данного пути
            for (int i = 0; i < this.Ntrc; ++i) {
                num_mesh1 = tracks[i].num_mesh;
                for (int j = 0; j < tracks[i].eps.length; ++j) {
                    if ((tracks[i].eps[j][0] > this.meshes[num_mesh1].X_end) || (tracks[i].eps[j][0] < this.meshes[num_mesh1].X_beg)) {
                        err_and_mes.data_error = true;
                        err_and_mes.messeg_data_error += "Путь номер "+(i+1)+": координаты точек ФОТ, ЗАЗ или сосред сопротивлен выходят за границы сетки";
                    }
                }
            }
        }
    }

    /** Метод для решения СЛАУ с 3диагонал ленточной матрицей
     методом двойной проходки
     на входе матрица и вектор правой части
     на выходе вектор ответов
     */
    public double[] solve_3diag_band(double[][] matrix_band, double[] vector_b) { // matrix_band – трёхдиагональная ленточная матрица, vector_b - вектор правой части
        int N = vector_b.length;
        double v[] = new double[N], u[] = new double[N], out[] = new double[N];
        //прямая проходка
        v[0] = matrix_band[2][1]/(-matrix_band[1][0]);
        u[0] = -vector_b[0]/(-matrix_band[1][0]);
        for (int i = 1; i < N-1; ++i) {
            v[i] = matrix_band[2][i+1]/(-matrix_band[1][i]-matrix_band[0][i-1]*v[i-1]);
            u[i] = (matrix_band[0][i-1]*u[i-1]-vector_b[i])/(-matrix_band[1][i]-matrix_band[0][i-1]*v[i-1]);
        }
        v[N-1] = 0;
        u[N-1] = (matrix_band[0][N-2]*u[N-2]-vector_b[N-1])/(-matrix_band[1][N-1]-matrix_band[0][N-2]*v[N-2]);
        //обратная проходка
        out[N-1] = u[N-1];
        for (int i = N-1; i > 0; --i) {
            out[i-1] = out[i]*v[i-1]+u[i-1];
        }
        return out;
    }

    /** Заполняет массивы номеров путей, сеток и индексов узлов для МПС
     в каждой строке два значения int: для начальной и конечной точки подключения
     */
    private void find_index_mps() {
        int N = this.mps.length;
        this.num_track_mps = new int[N][2]; // инициализируем массивы
        this.num_mesh_mps = new int[N][2];
        this.index_mps = new int[N][2];
        for (int i = 0; i < N; ++i) { //заполняем массивы
            this.num_track_mps[i][0] = (int) this.mps[i][0];
            this.num_track_mps[i][1] = (int) this.mps[i][1];
            this.num_mesh_mps[i][0] = this.tracks[num_track_mps[i][0]].num_mesh;
            this.num_mesh_mps[i][1] = this.tracks[num_track_mps[i][1]].num_mesh;
            this.index_mps[i][0] = this.meshes[this.num_mesh_mps[i][0]].find_near_index_over_mesh(this.mps[i][2]);
            this.index_mps[i][1] = this.meshes[this.num_mesh_mps[i][1]].find_near_index_over_mesh(this.mps[i][3]);
        }
    }

    /** Метод инициализирует ОТС
     */
    public void init_OTS() {
        int num_mesh; // номер сетки дял текущего пути

        verify_data(true, true, true); // проверка исходных данных полная
        for (int i = 0; i < this.Nmsh; ++i) {
            this.meshes[i].create_mesh(); //создаём сетку
        }
        find_index_mps(); // инициализируем массивы номеров начальной и конечной точки подключения для  МПС
        this.U_const = new double[this.mps.length]; // инициализируем массивы постоянной составляющей напряжения

        for (int i = 0; i < this.Ntrc; ++i) { // для обоих путей
            num_mesh = this.tracks[i].num_mesh;
            this.tracks[i].vector_b = new double[this.meshes[num_mesh].mesh_N]; // инициализируем вектор правой части
            this.meshes[num_mesh].create_3diag_matrix_band(tracks[i]); // рассчитываем 3диагональную ленточную матрицу
            this.tracks[i].u = new double[this.meshes[num_mesh].mesh_N]; // инициализация векторов напряжений в рельсах
            this.tracks[i].i = new double[this.meshes[num_mesh].mesh_N]; // токов в рельсах
            this.tracks[i].i_grnd = new double[this.meshes[num_mesh].mesh_N]; // токов из рельсов в землю
        }
        eval_a_x_find(); // рассчитываем матрицу коэффициентов влияния в МПС
    }

    /* Расчет постоянной составляющей напряжения в точках МПС
     */
    private void eval_U_const() {
        // Рассчитаем напряжение в узлах каждого пути от постоянных источников
        for (int i = 0; i < this.Ntrc; ++i) {
            this.tracks[i].vector_b = addValues_vector_b_1node(this.tracks[i].num_mesh, tracks[i].vector_b, tracks[i].fot, -1); // в точках ФОТ ток в одном узле с минусом
            this.tracks[i].vector_b = addValues_vector_b_2node(this.tracks[i].num_mesh, tracks[i].vector_b, tracks[i].eps, 1); // в точках ЭПС ток  в двух ближайших узлах
            this.tracks[i].u = solve_3diag_band(tracks[i].m3db, tracks[i].vector_b); //потенциал в рельсах в узлах сетки
        }

        // по параметрам МПС (номерам путей и координатам точек начала и конца) запишем в U1_const и U2_const  в узлах с МПС  напряжения в рельсах от постоянных источников
        for (int j = 0; j < this.mps.length; ++j) {
            this.U_const[j] = this.tracks[this.num_track_mps[j][0]].u[this.index_mps[j][0]]-this.tracks[this.num_track_mps[j][1]].u[this.index_mps[j][1]];
        }
    }

    /* Расчет коэффициентов влияния в МПС от их самих
     */
    private void eval_a_x_find() {
        int N_find = this.mps.length; // количестов поисковых точек
        double Imps = 1000, // уловный ток МПС для определения коэффициентов влияния
                u1[], u2[], // массивы напряжений в пути начальной и конечной точки подключения создаваемые током МПС
                a1_x_find[][] = new double[N_find][N_find], // временные массивы влияния в МПС в начальной и
                a2_x_find[][] = new double[N_find][N_find]; // в конечной
        this.a_x_find = new double[N_find][N_find]; // инициализируем общий массив влияния
        // Перебераем поисковые точки МПС
        for (int i = 0; i < N_find; ++i) {
            this.tracks[this.num_track_mps[i][0]].vector_b[this.index_mps[i][0]] = -Imps; // задаём ток в начальной точке подключения в данном пути
            this.tracks[this.num_track_mps[i][1]].vector_b[this.index_mps[i][1]] = Imps; // задаём ток в конечной точке подключения в данном пути
            u1 = solve_3diag_band(this.tracks[this.num_track_mps[i][0]].m3db, this.tracks[this.num_track_mps[i][0]].vector_b); //снимаем напряжение по сетке в пути для начальной точки подключения
            u2 = solve_3diag_band(this.tracks[this.num_track_mps[i][1]].m3db, this.tracks[this.num_track_mps[i][1]].vector_b); //снимаем напряжение по сетке в пути для конечной точки подключения
            tracks[this.num_track_mps[i][0]].vector_b[index_mps[i][0]] = 0; // обнуляем ток в текущей начальной точке подключения
            tracks[this.num_track_mps[i][1]].vector_b[index_mps[i][1]] = 0; // обнуляем ток в текущей конечной точке подключения
            for (int j = 0; j < N_find; ++j) { //снова проходим по всем точкам МПС
                if (this.num_track_mps[i][0] == this.num_track_mps[j][0]) { // если номер пути начальной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u1[this.index_mps[j][0]]/Imps; // находим коэффициенты влияния для текущей точки
                }
                if (this.num_track_mps[i][1] == this.num_track_mps[j][0]) { // если номер пути конечной точки МПС в цикле i =  номеру пути начальной точки МПС в цикле j
                    a1_x_find[i][j] = u2[this.index_mps[j][0]]/Imps; // находим коэффициенты влияния для текущей точки
                }
                if (this.num_track_mps[i][0] == this.num_track_mps[j][1]) { // если номер пути начальной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u1[this.index_mps[j][1]]/Imps; // находим коэффициенты влияния для текущей точки
                }
                if (this.num_track_mps[i][1] == this.num_track_mps[j][1]) { // если номер пути конечной точки МПС в цикле i =  номеру пути конечной точки МПС в цикле j
                    a2_x_find[i][j] = u2[this.index_mps[j][1]]/Imps; // находим коэффициенты влияния для текущей точки
                }
                this.a_x_find[i][j] = a1_x_find[i][j]-a2_x_find[i][j];
            }
        }
    }

    /** Заполняет вектор правой части по заданному двумерному массиву с коэффициентом тока
     для точек с одним узлом
     */
    private double[] addValues_vector_b_1node(int num_mesh, double vector_b[], double arr_XI[][], double coeff_I) { //num_mesh -номер сетки, vector_b -вектор правой части, coeff_I - коэффициент на который умножается значение тока
        for (int i = 0; i < arr_XI.length; ++i) {
            vector_b[this.meshes[num_mesh].find_near_index_over_mesh(arr_XI[i][0])] += arr_XI[i][1]*coeff_I;
        }
        return vector_b;
    }

    /** Добавляет значения в вектор правой части по заданному двумерному массиву с коэффициентом тока
     для точек с двумя узлами
     */
    private double[] addValues_vector_b_2node(int num_mesh, double vector_b[], double arr_XI[][], double coeff_I) { //num_mesh -номер сетки, vector_b - вектор правой части А; arr_XI - двумерный массив в каждой строке  два значения: 1 - координата км, 2 - ток А; coeff_I - коэффициент тока
        int index2[] = new int[2];
        double k, I1, I2;
        for (int i = 0; i < arr_XI.length; ++i) {
            index2 = this.meshes[num_mesh].find_2near_index_over_mesh(arr_XI[i][0]); // номер левого и правого узла
            k = (arr_XI[i][0]-this.meshes[num_mesh].X[0])/this.meshes[num_mesh].dX; // дробный индекс точки отнсительно номеров узлов сетки
            I2 = arr_XI[i][1]*coeff_I*(k-index2[0]); // ток в левый узел
            I1 = arr_XI[i][1]-I2; // ток в правый узел
            vector_b[index2[0]] += I1;
            vector_b[index2[1]] += I2;
        }
        return vector_b;
    }

    /** Заполняет вектор правой части всех точек нулями
     */
    private void zeros_vector_b() {
        int double_index[], //двойной индекс
                num_mesh; // номер сетки

        // проход по всем путям обнуление правой части ФОТ и ЭПС
        for (int j = 0; j < this.Ntrc; ++j) {
            // обнуление для точек в один узел сетки ФОТ
            num_mesh = this.tracks[j].num_mesh;
            for (int i = 0; i < this.tracks[j].fot.length; ++i) {
                this.tracks[j].vector_b[this.meshes[num_mesh].find_near_index_over_mesh(this.tracks[j].fot[i][0])] = 0;
            }
            //обнуление для точек в два узла сетки ЭПС
            for (int i = 0; i < this.tracks[j].eps.length; ++i) {
                double_index = this.meshes[num_mesh].find_2near_index_over_mesh(tracks[j].eps[i][0]);
                this.tracks[j].vector_b[double_index[0]] = 0;
                this.tracks[j].vector_b[double_index[1]] = 0;
            }
        }

        // проход по всем МПС
        for (int i = 0; i < this.mps.length; ++i) {
            this.tracks[this.num_track_mps[i][0]].vector_b[this.index_mps[i][0]] = 0;
            this.tracks[this.num_track_mps[i][1]].vector_b[this.index_mps[i][1]] = 0;
        }
    }

    /** Метод для расчёта мгновенной схемы ОТС
     с указанием начального значения тока в МПС
     */
    public boolean calc_ots(double init_I_poisk[]) {
        if (!verifi_I_poisk(init_I_poisk)) {
            return false;
        }
        return calc_I_poisk(init_I_poisk);
    }

    /** Метод для расчёта мгновенной схемы ОТС
     без указания начального значения тока в поисковых точках
     */
    public boolean calc_ots() {
        double init_I_poisk[] = new double[this.mps.length];
        return calc_I_poisk(init_I_poisk);
    }


    /**По сути рассчитывает токи в поисковых точках: МПС. После этого все граничные условия во всех точках с втекающим током определены
     Расчёт каждого поискового тока сводится к вычислению невязки напряжения на данном элементе и по величине невязки корректируется ток элемента
     По сути рассчитывается алгебраическая система уравнений методом Ньютона в цикле итераций
     */
    private boolean calc_I_poisk(double init_I_poisk[]) { // возвращает истина если расчёт сошёлся, ложь в обратном случае
        this.err_and_mes.reset_solver_error(); //обнулим ошибки решателя
        this.err_and_mes.calc_completed = true;
        verify_data(false, false, true); // проверка исходных данных только координ ЭПС, т.к. остальное проверено при инициализации
        if (this.err_and_mes.data_error) { //проверка если данные корректны
            this.err_and_mes.solver_error = true;
            this.err_and_mes.messeg_solver_error = "Расчёт невозможен. Ошибка исходных данных";
            return false;
        }

        int N_mps = this.mps.length, //количество МПС
                N_find = N_mps; //количество поисковых  всего
        double U_find[] = new double[N_find], // массивы напряжений на МПС, которые определяются всеми токами (известными и неизвестными)
                resid_U_mps[] = new double[N_mps], I_mps[] = init_I_poisk; // массив невязки напряжений на МПС и токов МПС
        double mean_resid, // средняя невязка
                limit_mean_resid = this.computing_settings.convergence_U, // задаём предельную невязку по достижении которой сходимость из класса computing_settings
                damping_factor = this.computing_settings.initial_damping_factor, //задаём коэффициент демпфирования текущее значение, на него умножается вычисленная по невязке напряжение корректирровка тока
                mean_resid_pred; // значение невязки по напряжению на предыдущем шаге итераций
        int iter = 0, iter_max = this.computing_settings.max_number_iterat; // счётчик итераций и максимальное число итераций
        boolean counter_not_exceeded, convergence_not_achieved; // непревышение итераций, недостижение сходимости - булевые переменные которые определяют выход из цикла итераций

        eval_U_const(); // рассчитываем напряжения на МПС от заданных токов ФОТ и ЭПС

        //нахождение токов в цикле итераций по невязке напряжения на МПС
        counter_not_exceeded = true;
        convergence_not_achieved = true;
        mean_resid = 100; //начальное значение средняя невязка до первой итерации
        while (counter_not_exceeded && convergence_not_achieved) {
            mean_resid_pred = mean_resid; //предыдущая невязка обновление
            mean_resid = 0; //текущая невязка скидывается
            for (int i = 0; i < N_find; ++i) {
                U_find[i] = U_const[i]; //начинаем с постоянного напряжения (от заданных источников тока ФОТ ЭПС)
                for (int j = 0; j < N_mps; ++j) {
                    U_find[i] += I_mps[j]*this.a_x_find[i][j]; //добавляем напряжение от МПС
                }
                resid_U_mps[i] = U_find[i]-I_mps[i]*this.mps[i][4]; //невязка напряжения на МПС, уже изветны напряжения и Р1 и Р2
                I_mps[i] += damping_factor*resid_U_mps[i]/(-0.5*this.a_x_find[i][i]+this.mps[i][4]); //корректируем текущий поисковый ток пропорционально невязке по напряжению в этом элементе с учётом коэф. демпфирования
                mean_resid += Math.abs(resid_U_mps[i]); //обновляем невязку
            }
            mean_resid = mean_resid/N_find; //невязка именно средняя
            //если после первой итерации возрастает средняя невязка mean_resid по сравнению с ней же на предыдущей итерации mean_resid_pred, то коэффициент демпфирования в методе Ньютона уменьшаем в 0.7 раз
            if (iter > 0) {
                if (mean_resid > mean_resid_pred) {
                    damping_factor = damping_factor*0.7;
                }
            }
            iter += 1; //обновляем счётчик итераций
            counter_not_exceeded = (iter < iter_max); // обновляем булевые переменные выхода из цикла итераций
            convergence_not_achieved = (mean_resid > limit_mean_resid);
            //debugLog("iter="+iter+" ,mean_resid="+mean_resid+", damping_factor="+damping_factor);
        }
        //debugLog("iter="+iter);
        this.computing_settings.current_state_solver = new double[]{iter-1, mean_resid, damping_factor}; // записываем текущее состояние решателя
        this.I_mps = I_mps; // заносим токи в МПС в массивы родительского класса
        eval_node_from_all_I();
        zeros_vector_b();

        if (convergence_not_achieved) {
            this.err_and_mes.solver_error = true;
            this.err_and_mes.messeg_solver_error = "Превышено максимальное число итераций равное "+this.computing_settings.max_number_iterat+" заданная сходимость (сред невязка по напряжен) равная "+this.computing_settings.convergence_U+" не достигнута.";
        }
        return !convergence_not_achieved; // возврат обратное к несходимости: истина если расчёт сошёлся
    }

    /*При найденных токах МПС рассчитывает напряжения и ток в рельсах решением М3ДЛ
     */
    private void eval_node_from_all_I() {
        double g_rh, g_lf; //условная проводимость слева и справа от узла сетки на схеме дискретизации рельсов
        int num_mesh, // номер сетки для начального и конечного пути МПС
                N; // количество узлов сетки

        // добавляем в вектор правой часть токи МПС (токи ФОТ и ЭПС добавлены процедурой eval_U_const() )
        for (int i = 0; i < this.mps.length; ++i) {
            this.tracks[this.num_track_mps[i][0]].vector_b[this.index_mps[i][0]] -= this.I_mps[i];
            this.tracks[this.num_track_mps[i][1]].vector_b[this.index_mps[i][1]] += this.I_mps[i];
        }

        //потенциал в рельсах в узлах сетки для каждого пути
        for (int i = 0; i < this.Ntrc; ++i) {
            this.tracks[i].u = solve_3diag_band(this.tracks[i].m3db, this.tracks[i].vector_b);
        }

        //расчёт токов в рельсах и тока в земле в узлах сетки для каждого пути
        for (int j = 0; j < this.Ntrc; ++j) {
            num_mesh = this.tracks[j].num_mesh; //номер сетки для текущего пути
            N = this.meshes[num_mesh].mesh_N; // число узлов сетки для текущего пути
            //ток в рельсах и земле для первого узла
            g_lf = 1/this.tracks[j].Rv0; //проводимость слева от первого узла через волновое сопротивление в начале
            g_rh = -this.tracks[j].m3db[0][0]; //проводимость справа от первого узла через нижнюю диагональ первый элемент
            this.tracks[j].i[0] = 0.5*((0-this.tracks[j].u[0])*g_lf+(this.tracks[j].u[0]-this.tracks[j].u[1])*g_rh); //ток первого узла как полусумма токов слева и справа от него
            this.tracks[j].i_grnd[0] = this.tracks[j].u[0]*(this.tracks[j].m3db[1][0]+this.tracks[j].m3db[0][0]); //ток в земле для первого узла
            //ток в рельсах и земле для остальных узлов со второго до предпоследнего
            for (int i = 1; i < N-1; ++i) {
                g_lf = -this.tracks[j].m3db[0][i-1]; //проводимость слева от узла через нижнюю диагональ
                g_rh = -this.tracks[j].m3db[2][i+1]; //проводимость справа от узла через верхнюю диагональ
                this.tracks[j].i[i] = 0.5*((this.tracks[j].u[i-1]-this.tracks[j].u[i])*g_lf+(this.tracks[j].u[i]-this.tracks[j].u[i+1])*g_rh); //ток  узла как полусумма токов слева и справа от него
                this.tracks[j].i_grnd[i] = this.tracks[j].i_grnd[i-1]+this.tracks[j].u[i]*(this.tracks[j].m3db[1][i]+this.tracks[j].m3db[0][i-1]+this.tracks[j].m3db[2][i+1]); //ток в земле для узла
            }
            //ток в рельсах для последнего узла
            g_lf = -this.tracks[j].m3db[2][N-1]; //проводимость слева от последнего узла через верхнюю диагональ последний элемент
            g_rh = 1/this.tracks[j].Rvn; //проводимость справа от последнего узла через волновое сопротивление в конце
            this.tracks[j].i[N-1] = 0.5*((this.tracks[j].u[N-2]-this.tracks[j].u[N-1])*g_lf+(this.tracks[j].u[N-1]-0)*g_rh); // ток в рельсе для псоледнего узла
            this.tracks[j].i_grnd[N-1] = this.tracks[j].i_grnd[N-2]+this.tracks[j].u[N-1]*(this.tracks[j].m3db[1][N-1]+this.tracks[j].m3db[2][N-1]-1/this.tracks[j].Rvn); //ток в земле для последнего узла
        }
    }


    //------------------------------------------------------Геттеры и сеттеры головного класса + внутренние процедуры для них-------------------------------------------------------------------------

    // возвращает параметр распределенный по узлам сетки в привязке к координатам
    private double[] return_X_param_rail(int num_mesh, double X_rail[], double[] arr_param_node) { //num_mesh - номер сетки, X_rail[] - массив координат, км; arr_param_node - массив параметра рсапредленный по узлам сетки
        int N = X_rail.length, indexes[];
        double out[] = new double[N],
                proportion_node_left;
        //проверки при некорректности возврат null
        if (this.err_and_mes.data_error || !this.err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            return null;
        }
        if ((X_rail[0] < this.meshes[num_mesh].X[0]) || (X_rail[N-1] > this.meshes[num_mesh].X[this.meshes[num_mesh].mesh_N-1])) { // проверка если координаты расчётного массива за пределами сетки
            return null;
        }
        //непосредственно распредление параметра в узлах сетки на массив координат в пределах сетки
        for (int i = 0; i < N; ++i) {
            indexes = this.meshes[num_mesh].find_2near_index_over_mesh(X_rail[i]); // индексы левого и правого узла сетки относительно координаты X_rail[i]
            proportion_node_left = (this.meshes[num_mesh].X[indexes[1]]-X_rail[i])/this.meshes[num_mesh].dX; // доля для левого узла сетки обратно пропорциональна расстоянию до правого узла
            out[i] = arr_param_node[indexes[0]]*proportion_node_left+arr_param_node[indexes[1]]*(1-proportion_node_left); // напряжение в левом узле на его долю + напряжение в правом узле на его долю
        }
        return out;
    }

    //геттер возвращает напряжение рельс-земля в виде массива по узлам сетки
    public double[] get_U_rail(int num_track) { // num_track - номер пути
        return this.tracks[num_track].u;
    }

    //геттер возвращает ток в рельсах в виде массива по узлам сетки
    public double[] get_I_rail(int num_track) { // num_track - номер пути
        return this.tracks[num_track].i;
    }

    //геттер возвращает ток в земле в виде массива по узлам сетки для номеров путей заданных массивом
    // все пути для которых вычисляется ток в земле должны быть с единой сеткой
    public double[] get_I_grnd(int nums_tracks[]) { // nums_tracks - массив содержащий номера путей для которых суммировать ток в земле
        int num_mesh = this.tracks[nums_tracks[0]].num_mesh, // номер сетки для первого пути в массиве
                N = this.meshes[num_mesh].mesh_N; // число узлов сетки
        double i_grnd_sum[] = new double[N];
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < nums_tracks.length; ++i) {
                i_grnd_sum[j] += this.tracks[nums_tracks[i]].i_grnd[j];
            }
        }
        return i_grnd_sum;
    }

    //геттер возвращает ток в земле в виде массива по узлам сетки в сумме для всех путей
    public double[] get_I_grnd() {
        int nums_tracks[] = new int[this.Ntrc];
        for (int i = 0; i < this.Ntrc; ++i) {
            nums_tracks[i] = i+1;
        }
        return get_I_grnd(nums_tracks);
    }

    //геттер возвращает напряжение рельс-земля в виде массива по заданным координатам точек
    public double[] get_U_rail(int num_track, double X_rail[]) { //  num_track - номер пути; X_rail - массив координат км;
        int num_mesh = this.tracks[num_track].num_mesh; // номер сетки для пути
        return return_X_param_rail(num_mesh, X_rail, get_U_rail(num_track));
    }

    //геттер возвращает ток в рельсах в виде массива по заданным координатам точек
    public double[] get_I_rail(int num_track, double X_rail[]) { // num_track - номер пути; X_rail - массив координат км;
        int num_mesh = this.tracks[num_track].num_mesh; // номер сетки для пути
        return return_X_param_rail(num_mesh, X_rail, get_I_rail(num_track));
    }

    //геттер возвращает ток в земле в виде массива для заданных путей по заданным координатам точек
    public double[] get_I_grnd(int nums_tracks[], double X_rail[]) { // nums_tracks - массив содержащий номера путей для которых суммировать ток в земле; X_rail - массив координат км
        int num_mesh = this.tracks[nums_tracks[0]].num_mesh; // номер сетки для первого пути в массиве
        return return_X_param_rail(num_mesh, X_rail, get_I_grnd(nums_tracks));
    }


    //возвращает массив с токами в МПС
    public double[] get_I_poisk() {
        if (this.err_and_mes.data_error || !this.err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим пустой массив
            return null;
        }
        return this.I_mps;
    }

    /*проверка коорректности входного двумерного массива I_poisk[] содержащего токи МПС
     возвращает true если корректен, и в противном случае false
     */
    private boolean verifi_I_poisk(double[] I_poisk) {
        if (I_poisk.length != this.mps.length) {
            this.err_and_mes.messeg_solver_error = "Ошибка при вводе массива токов МПС. Количество элементов массива должно совпадать с количеством соответствующих точек МПС";
            return false;
        }
        return true;
    }

    /*сеттер задает значения в МПС без выполнения непосредственного системы уравнений
      возвращает true при успешном выполнении процедуры, false  - в противном случае
      на вход принимает двумерный массив содержащий три строки, это массивы токов в следующем порядке: МПС, ЗАЗ1, ЗАЗ2
    */
    public boolean set_I_poisk_no_call(double[] I_poisk) {
        if (this.err_and_mes.data_error) { // проверка если ошибка исходных данных, то выводим  false. Без кооректного задания ФОТ, ЭПС1 и ЭПС2 - расчёт даже в этом случае невозможен
            this.err_and_mes.messeg_solver_error = "Задание тока в МПС не возможно. Ошибка исходных данных";
            return false;
        }
        if (!verifi_I_poisk(I_poisk)) { //вызываем метод проверки корректности введенного двумерного массива I_poisk
            return false;
        }
        //записываем массивы токов в поисковых точках из заданных на входе
        this.I_mps = I_poisk;
        zeros_vector_b(); // обнуление вектора правой части всех путей
        eval_U_const(); // расчёт от постоянных источников
        eval_node_from_all_I(); // расчёт всех узлов с учётом токов МПС
        return true;
    }

    //расчитывает и выдаёт мгновенную мощность потерь в ОТС (все элементы ОТС) в кВт
    public double get_P_ots() {
        if (this.err_and_mes.data_error || !this.err_and_mes.calc_completed) { // проверка если ошибка исходных данных или расчёт не выполнен, то выводим  -1
            return -1;
        }
        double P_ots = 0; // мощность ОТС накапливается в этой переменной
        for (int i = 0; i < this.Ntrc; ++i) { // проход по всем путям
            for (int j = 0; j < this.tracks[i].fot.length; ++j) { // проход по всем ФОТ данного пути, его ток с минусов
                P_ots += -this.tracks[i].fot[j][1]*get_U_rail(i, new double[]{this.tracks[i].fot[j][0]})[0];
            }
            for (int j = 0; j < this.tracks[i].eps.length; ++j) { // проход по всем ЭПС данного пути
                P_ots += this.tracks[i].eps[j][1]*get_U_rail(i, new double[]{this.tracks[i].eps[j][0]})[0];
            }
        }
        return P_ots;
    }
}