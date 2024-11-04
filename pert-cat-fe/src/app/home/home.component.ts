import { Component, OnInit } from '@angular/core';
import {HeaderComponent} from "../shared/header/header.component";
import {MatCardModule} from "@angular/material/card";
import {FooterComponent} from "../shared/footer/footer.component";
import { BigQueryService } from '../bigquery.service';
import { MatTableModule } from '@angular/material/table';
import { CommonModule } from '@angular/common';


@Component({
  selector: 'app-home',
  standalone: true,
  imports: [HeaderComponent, FooterComponent, MatCardModule, MatTableModule, CommonModule],
  templateUrl: './home.component.html',
  styleUrl: './home.component.css'
})
export class HomeComponent implements OnInit {
  data: any[] = [];

  constructor(private bigQueryService: BigQueryService) { }

  ngOnInit() {
    this.fetchData();
  }

  fetchData() {
    this.bigQueryService.getData().subscribe({
      next: (response) => {
        this.data = response.rows;
      },
      error: (error) => {
        console.error('Error fetching data:', error);
      }
    });
  }
}
