import { Routes } from '@angular/router';
import {HomeComponent} from "./home/home.component";
import {DataPortalComponent} from "./data-portal/data-portal.component";
import {DetailsComponent} from "./data-portal/details/details.component";
import {DashboardsComponent} from "./dashboards/dashboards.component";
import {DataAnalyticsComponent} from "./data-analytics/data-analytics.component";
import {AboutComponent} from "./about/about.component";

export const routes: Routes = [
  { path: 'home', title: 'Perturbation Catalogue | Home', component: HomeComponent },
  { path: '', redirectTo: 'home', pathMatch: 'full' },
  { path: 'data_portal', title: 'Perturbation Catalogue | Data Portal', component: DataPortalComponent },
  { path: 'data_portal/:urn', title: 'Perturbation Catalogue | Details', component: DetailsComponent },
  { path: 'dashboards', title: 'Perturbation Catalogue | Dashboards', component: DashboardsComponent },
  { path: 'data-analytics', title: 'Perturbation Catalogue | Data Analytics', component: DataAnalyticsComponent },
  { path: 'about', title: 'Perturbation Catalogue | About', component: AboutComponent },
];
