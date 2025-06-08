#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

class PopulationGrowthModel {
private:
    double r;        // Laju pertumbuhan intrinsik
    double K;        // Carrying capacity 
    double dt;       // step size
    double t_max;    // time maksimum simulasi
    
public:
    PopulationGrowthModel(double growth_rate, double carrying_capacity, double step_size, double max_time) 
        : r(growth_rate), K(carrying_capacity), dt(step_size), t_max(max_time) {}
    
    // Model Pertumbuhan Eksponensial dengan dP/dt = r*P
    double exponential_growth(double P) {
        return r * P;
    }
    
    // Model Pertumbuhan Logistik yang dngan dP/dt = r*P*(1 - P/K)
    double logistic_growth(double P) {
        return r * P * (1.0 - P/K);
    }
    
    // Metode Euler untuk menyelesaikan ODE
    vector<pair<double, double>> solve_euler(double P0, string model_type) {
        vector<pair<double, double>> results;
        double t = 0.0;
        double P = P0;
        
        while (t <= t_max) {
            results.push_back({t, P});
            
            double dPdt;
            if (model_type == "exponential") {
                dPdt = exponential_growth(P);
            } else {
                dPdt = logistic_growth(P);
            }
            
            P = P + dt * dPdt;
            t = t + dt;
        }
        
        return results;
    }
    
    // Metode Runge-Kutta Order 4 (RK4)
    vector<pair<double, double>> solve_rk4(double P0, string model_type) {
        vector<pair<double, double>> results;
        double t = 0.0;
        double P = P0;
        
        while (t <= t_max) {
            results.push_back({t, P});
            
            double k1, k2, k3, k4;
            
            if (model_type == "exponential") {
                k1 = dt * exponential_growth(P);
                k2 = dt * exponential_growth(P + k1/2);
                k3 = dt * exponential_growth(P + k2/2);
                k4 = dt * exponential_growth(P + k3);
            } else {
                k1 = dt * logistic_growth(P);
                k2 = dt * logistic_growth(P + k1/2);
                k3 = dt * logistic_growth(P + k2/2);
                k4 = dt * logistic_growth(P + k3);
            }
            
            P = P + (k1 + 2*k2 + 2*k3 + k4) / 6.0;
            t = t + dt;
        }
        
        return results;
    }
    
    // Metode Prediksi-Koreksi (Adams-Bashforth-Moulton)
    vector<pair<double, double>> solve_predictor_corrector(double P0, string model_type) {
        vector<pair<double, double>> results;
        double t = 0.0;
        double P = P0;
        
        // Gunakan RK4 untk langkah pertama
        results.push_back({t, P});
        
        double dPdt_prev;
        if (model_type == "exponential") {
            dPdt_prev = exponential_growth(P);
        } else {
            dPdt_prev = logistic_growth(P);
        }
        
        P = P + dt * dPdt_prev;
        t = t + dt;
        results.push_back({t, P});
        
        // Lanjutkan dengan metode prediksi-koreksi
        while (t <= t_max) {
            double dPdt_curr;
            if (model_type == "exponential") {
                dPdt_curr = exponential_growth(P);
            } else {
                dPdt_curr = logistic_growth(P);
            }
            
            // Prediksi (Adams-Bashforth)
            double P_pred = P + dt * (1.5 * dPdt_curr - 0.5 * dPdt_prev);
            
            // Koreksi (Adams-Moulton)
            double dPdt_pred;
            if (model_type == "exponential") {
                dPdt_pred = exponential_growth(P_pred);
            } else {
                dPdt_pred = logistic_growth(P_pred);
            }
            
            P = P + dt * (0.5 * dPdt_pred + 0.5 * dPdt_curr);
            t = t + dt;
            
            results.push_back({t, P});
            dPdt_prev = dPdt_curr;
        }
        
        return results;
    }
    
    // Fungsi untuk menghitung solusi analitik yaitu jika tersedia
    double analytical_exponential(double P0, double t) {
        return P0 * exp(r * t);
    }
    
    double analytical_logistic(double P0, double t) {
        return (K * P0) / (P0 + (K - P0) * exp(-r * t));
    }
    
    // Fungsi untuk menghitung error
    double calculate_rmse(vector<pair<double, double>>& numerical, vector<pair<double, double>>& analytical) {
        double sum_squared_error = 0.0;
        int n = min(numerical.size(), analytical.size());
        
        for (int i = 0; i < n; i++) {
            double error = numerical[i].second - analytical[i].second;
            sum_squared_error += error * error;
        }
        
        return sqrt(sum_squared_error / n);
    }
    
    // Fungsi untuk menyimpan hasil ke file CSV
    void save_to_csv(const vector<pair<double, double>>& data, const string& filename) {
        ofstream file(filename);
        file << "Time,Population\n";
        
        for (const auto& point : data) {
            file << fixed << setprecision(4) << point.first << "," << point.second << "\n";
        }
        
        file.close();
        cout << "Data saved to " << filename << endl;
    }
};

// Fungsi untuk menampilkan menu
void display_menu() {
    cout << "\n=== APLIKASI PEMODELAN PERTUMBUHAN POPULASI ===" << endl;
    cout << "1. Model Pertumbuhan Eksponensial" << endl;
    cout << "2. Model Pertumbuhan Logistik" << endl;
    cout << "3. Perbandingan Metode Numerik" << endl;
    cout << "4. Keluar" << endl;
    cout << "Pilih opsi (1-4): ";
}

// Fungsi untuk menampilkan hasil
void display_results(const vector<pair<double, double>>& results, const string& method_name) {
    cout << "\n--- Hasil " << method_name << " ---" << endl;
    cout << setw(10) << "Waktu" << setw(15) << "Populasi" << endl;
    cout << string(25, '-') << endl;
    
    for (size_t i = 0; i < results.size(); i += max(1, (int)results.size()/10)) {
        cout << setw(10) << fixed << setprecision(2) << results[i].first 
             << setw(15) << fixed << setprecision(0) << results[i].second << endl;
    }
}

int main() {
    cout << fixed << setprecision(2);
    
    int choice;
    do {
        display_menu();
        cin >> choice;
        
        switch (choice) {
            case 1: {
                // Model Pertumbuhan Eksponensial
                cout << "\n=== MODEL PERTUMBUHAN EKSPONENSIAL ===" << endl;
                cout << "Persamaan: dP/dt = r*P" << endl;
                
                double r, P0, dt, t_max;
                cout << "Masukkan laju pertumbuhan (r): ";
                cin >> r;
                cout << "Masukkan populasi awal (P0): ";
                cin >> P0;
                cout << "Masukkan step size (dt): ";
                cin >> dt;
                cout << "Masukkan waktu maksimum: ";
                cin >> t_max;
                
                PopulationGrowthModel model(r, 0, dt, t_max);
                
                // Hitung dengan berbagai metode
                auto euler_results = model.solve_euler(P0, "exponential");
                auto rk4_results = model.solve_rk4(P0, "exponential");
                
                // Hitung solusi analitik
                vector<pair<double, double>> analytical_results;
                for (double t = 0; t <= t_max; t += dt) {
                    analytical_results.push_back({t, model.analytical_exponential(P0, t)});
                }
                
                display_results(euler_results, "Metode Euler");
                display_results(rk4_results, "Metode RK4");
                display_results(analytical_results, "Solusi Analitik");
                
                // Hitung error
                double rmse_euler = model.calculate_rmse(euler_results, analytical_results);
                double rmse_rk4 = model.calculate_rmse(rk4_results, analytical_results);
                
                cout << "\n--- Analisis Error ---" << endl;
                cout << "RMSE Euler: " << scientific << rmse_euler << endl;
                cout << "RMSE RK4: " << scientific << rmse_rk4 << endl;
                
                // Simpan ke file
                model.save_to_csv(euler_results, "exponential_euler.csv");
                model.save_to_csv(rk4_results, "exponential_rk4.csv");
                model.save_to_csv(analytical_results, "exponential_analytical.csv");
                
                break;
            }
            
            case 2: {
                // Model Pertumbuhan Logistik
                cout << "\n=== MODEL PERTUMBUHAN LOGISTIK ===" << endl;
                cout << "Persamaan: dP/dt = r*P*(1 - P/K)" << endl;
                
                double r, K, P0, dt, t_max;
                cout << "Masukkan laju pertumbuhan (r): ";
                cin >> r;
                cout << "Masukkan carrying capacity (K): ";
                cin >> K;
                cout << "Masukkan populasi awal (P0): ";
                cin >> P0;
                cout << "Masukkan step size (dt): ";
                cin >> dt;
                cout << "Masukkan waktu maksimum: ";
                cin >> t_max;
                
                PopulationGrowthModel model(r, K, dt, t_max);
                
                // Hitung dengan berbagai metode
                auto euler_results = model.solve_euler(P0, "logistic");
                auto rk4_results = model.solve_rk4(P0, "logistic");
                auto pc_results = model.solve_predictor_corrector(P0, "logistic");
                
                // Hitung solusi analitik
                vector<pair<double, double>> analytical_results;
                for (double t = 0; t <= t_max; t += dt) {
                    analytical_results.push_back({t, model.analytical_logistic(P0, t)});
                }
                
                display_results(euler_results, "Metode Euler");
                display_results(rk4_results, "Metode RK4");
                display_results(pc_results, "Metode Prediksi-Koreksi");
                display_results(analytical_results, "Solusi Analitik");
                
                // Hitung error
                double rmse_euler = model.calculate_rmse(euler_results, analytical_results);
                double rmse_rk4 = model.calculate_rmse(rk4_results, analytical_results);
                double rmse_pc = model.calculate_rmse(pc_results, analytical_results);
                
                cout << "\n--- Analisis Error ---" << endl;
                cout << "RMSE Euler: " << scientific << rmse_euler << endl;
                cout << "RMSE RK4: " << scientific << rmse_rk4 << endl;
                cout << "RMSE Prediksi-Koreksi: " << scientific << rmse_pc << endl;
                
                // Simpan ke file
                model.save_to_csv(euler_results, "logistic_euler.csv");
                model.save_to_csv(rk4_results, "logistic_rk4.csv");
                model.save_to_csv(pc_results, "logistic_predictor_corrector.csv");
                model.save_to_csv(analytical_results, "logistic_analytical.csv");
                
                break;
            }
            
            case 3: {
                // Perbandingan Metode Numerik
                cout << "\n=== PERBANDINGAN METODE NUMERIK ===" << endl;
                cout << "Menggunakan data default untuk perbandingan" << endl;
                
                // Parameter default
                double r = 0.1, K = 1000, P0 = 50, t_max = 50;
                vector<double> step_sizes = {0.1, 0.5, 1.0, 2.0};
                
                cout << setw(10) << "Step Size" << setw(15) << "RMSE Euler" 
                     << setw(15) << "RMSE RK4" << setw(20) << "RMSE Pred-Corr" << endl;
                cout << string(60, '-') << endl;
                
                for (double dt : step_sizes) {
                    PopulationGrowthModel model(r, K, dt, t_max);
                    
                    auto euler_results = model.solve_euler(P0, "logistic");
                    auto rk4_results = model.solve_rk4(P0, "logistic");
                    auto pc_results = model.solve_predictor_corrector(P0, "logistic");
                    
                    vector<pair<double, double>> analytical_results;
                    for (double t = 0; t <= t_max; t += dt) {
                        analytical_results.push_back({t, model.analytical_logistic(P0, t)});
                    }
                    
                    double rmse_euler = model.calculate_rmse(euler_results, analytical_results);
                    double rmse_rk4 = model.calculate_rmse(rk4_results, analytical_results);
                    double rmse_pc = model.calculate_rmse(pc_results, analytical_results);
                    
                    cout << setw(10) << dt << setw(15) << scientific << rmse_euler 
                         << setw(15) << scientific << rmse_rk4 
                         << setw(20) << scientific << rmse_pc << endl;
                }
                
                cout << "\nKesimpulan:" << endl;
                cout << " Metode RK4 umumnya memberikan akurasi terbaik" << endl;
                cout << "Error meningkat seiring dengan bertambahnya step size" << endl;
                cout << "Metode Prediksi-Koreksi memberikan hasil yang baik dengan komputasi efisien" << endl;
                
                break;
            }
            
            case 4:
                cout << "Terima kasih telah gunain aplikasi!" << endl;
                break;
                
            default:
                cout << "Pilihan tidak valid!" << endl;
                break;
        }
        
    } while (choice != 4);
    
    return 0;
}