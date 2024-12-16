import { Routes } from '@angular/router';
import {HomeComponent} from "./home/home.component";
import {DashboardComponent} from "./dashboard/dashboard.component";
import {SearchComponent} from "./search/search.component";
import {DataPortalComponent} from "./data-portal/data-portal.component";

export const routes: Routes = [
  { path: 'home', title: 'Perturbation Catalogue | Home',  component: HomeComponent },
  { path: '', redirectTo: 'home', pathMatch: 'full' },
  { path: 'data_portal', title: 'Perturbation Catalogue | Data Portal', component: DataPortalComponent },
  { path: 'dashboards', title: 'Perturbation Catalogue | Dashboards', component: DashboardComponent },
  { path: 'search', title: 'Perturbation Catalogue | Search', component: SearchComponent },
];
