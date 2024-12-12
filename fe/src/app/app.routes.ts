import { Routes } from '@angular/router';
import {HomeComponent} from "./home/home.component";
import {DashboardComponent} from "./dashboard/dashboard.component";

export const routes: Routes = [
  { path: 'home', title: 'Perturbation Catalogue | Home',  component: HomeComponent },
  { path: '', redirectTo: 'home', pathMatch: 'full' },
  { path: 'dashboards', component: DashboardComponent },
];
