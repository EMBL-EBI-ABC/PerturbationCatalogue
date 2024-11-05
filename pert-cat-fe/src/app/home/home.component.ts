import { Component, OnInit } from '@angular/core';
import { HeaderComponent } from "../shared/header/header.component";
import { MatCardModule } from "@angular/material/card";
import { FooterComponent } from "../shared/footer/footer.component";
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
  data: any;

  constructor(
    private elasticService: ElasticService,
  ) { };


  ngOnInit() {
    this.fetchData();
  }

  getKeys(obj: any): string[] {
    return Object.keys(obj);
  }

  fetchData() {
    this.elasticService.getData().subscribe(
      (response) => {
        this.data = response;
        console.log('Data received:', this.data.length);
      },
      (error) => {
        console.error('Error fetching data:', error);
      }
    );
  }
}
