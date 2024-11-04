import { Component, OnInit } from '@angular/core';
import {HeaderComponent} from "../shared/header/header.component";
import {MatCardModule} from "@angular/material/card";
import {FooterComponent} from "../shared/footer/footer.component";
import { ElasticService } from '../elastic.service';
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

  constructor(private elasticService: ElasticService) { }

  ngOnInit() {
    this.fetchData();
  }

  fetchData() {
    this.elasticService.getData().subscribe({
      next: (response) => {
        this.data = response.rows;
      },
      error: (error) => {
        console.error('Error fetching data:', error);
      }
    });
  }
}
